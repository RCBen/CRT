function [rho, rhoH2, rhoId, rhoG] = getRho(const, pr, Tr, x18, x0)
    % [rho, rhoH2, rhoId, rhoG] = getRho(const, pr, Tr, x18, [x0])
    %
    % density of liquid mixture, containing H2 (g l^-1)
    %            liquid H2                     (g l^-1)
    %            ideal gas mixture             (g l^-1)
    %            real gas mixture              (g l^-1)
    % 
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % x18 / x0 (-)      scalar / vector / matrix (nargout < 4)
    %
    % see DOI: 10.1002/aic.690190120
    %    ISBN: 0-07-0011682-2 (pp. 5.5 ff & 5.23 ff & 6.30 f)
    
    if nargin <= 4
        x0 = 1 - x18;                                                       % molar fraction of H0-DBT in liquid phase (-)
    else
        if ~isequal(size(x18), size(x0))
            error('size(x18) ~= size(x0)')
        end
    end
    xH2 = 1 - x18 - x0;                                                     % molar fraction of H2 in liquid phase (-)
    HG = x18./(1-xH2);                                                      % hydrogenation grade (-)
    
    nP = numel(pr);
    nTheta = numel(Tr);
    if ~isempty(pr)
        if ~isrow(pr)
            error('pr has to be a scalar or a row-vector')
        elseif isrow(pr) && (nP ~= nTheta || ~isscalar(pr) && pr(1) ~= pr(2))
            pr = repelem(pr,1,nTheta);
            Tr = repmat(Tr,1,nP);
        end
    end
    if ~isrow(Tr) && ~isequal(size(HG), size(Tr))
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
            x18 = x18';
            x0 = x0';
            xH2 = xH2';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    elseif nargout == 4 && ~isvector(HG)
        error('HG has to be a vector') 
    end
          
    %% experimental data (t = 20, 30, 40, 50, 60, 70 °C)   
    rho = (1.2516 - 7.1058e-04*Tr - 1.5056e-01*HG + ...                     % density of liquid mixture, wo H2 (g l^-1)
        7.4948e-05*Tr.*HG)*1000;
    
    %% sum of liquid partial molar volumes of H2 & H2-free liquid solvent
    % see DOI: 10.1002/aic.690190120
    %    ISBN: 0-07-0011682-2 (pp. 5.5 ff & 5.23 ff & 6.30 f)
    if nargout >= 2 || any(xH2(:) ~= 0)
        [gamma18, gamma0] = RST(const, Tr, x18, x0);                        % activity coefficient of H18-DBT & H0-DBT (-)
        [pLV18, pLV0] = getPLV(const, Tr);                                  % vapour pressure of H18-DBT & H0-DBT (bar)  
        
        M = HG*const.M18 + (1-HG)*const.M0;                                 % molar mass of liquid mixture, wo H2 (g mol^-1)        
        v = M./rho;                                                         % molar volume of liquid mixture, wo H2 (l mol^-1)
        pLV = HG.*gamma18.*pLV18 + (1-HG).*gamma0.*pLV0;                    % vapour pressure of liquid mixture, wo H2 (bar)
        
        % mixing rules of Chueh & Prausnitz
        k1818 = 1 - (sqrt(const.vc18^(1/3)*const.vc18^(1/3))/ ...
            ((const.vc18^(1/3) + const.vc18^(1/3))/2))^3;
        k180 = 1 - (sqrt(const.vc18^(1/3)*const.vc0^(1/3))/ ...
            ((const.vc18^(1/3) + const.vc0^(1/3))/2))^3;	
        k018 = 1 - (sqrt(const.vc0^(1/3)*const.vc18^(1/3))/ ...
            ((const.vc0^(1/3) + const.vc18^(1/3))/2))^3;
        k00 = 1 - (sqrt(const.vc0^(1/3)*const.vc0^(1/3))/ ...
            ((const.vc0^(1/3) + const.vc0^(1/3))/2))^3;

        Tc1818 = sqrt(const.Tc18*const.Tc18)*(1-k1818);
        Tc180 = sqrt(const.Tc18*const.Tc0)*(1-k180);
        Tc018 = sqrt(const.Tc0*const.Tc18)*(1-k018);
        Tc00 = sqrt(const.Tc0*const.Tc0)*(1-k00);

        phi18 = HG*const.vc18^(2/3)./ ...                                   % surface fraction of H18-DBT (-)
            (HG*const.vc18^(2/3) + (1-HG)*const.vc0^(2/3));           
        phi0 = 1 - phi18;                                                   % surface fraction of H0-DBT (-)       

        Tc = phi18.^2*Tc1818 + phi18.*phi0*Tc180 + ...                      % critical temperature of liquid mixture, wo H2 (K)
            phi0.*phi18*Tc018 + phi0.^2*Tc00;
        vc = (HG*const.vc18^(2/3) + (1-HG)*const.vc0^(2/3)).^(3/2);         % critical molar volume of liquid mixture, wo H2 (l mol^-1)
        Zc = HG*const.Zc18 + (1-HG)*const.Zc0;                              % critical compressibility factor of liquid mixture , wo H2 (-)
        pc = (Zc*const.R.*Tc./(vc/1000))*1e-05;                             % critical pressure of liquid mixture, wo H2 (bar)

        w = HG*const.w18 + (1-HG)*const.w0;                                 % acentric factor of liquid mixture, wo H2 (-)  
        
        A = -0.9556 + 1.244*w;
        B = 1.712 + 0.283*w;
        C = 1.282 + 0.840*w;
        D = -4.40 + 3.10*sqrt(w);
        E = 4.20 - 3.10*sqrt(w);

        vH2 = 1./ ...                                                       % molar volume of H2 in liquid phase (l mol^-1)
            (1./(const.vcH2*exp(A + B.*(Tr./Tc).^3 + C.*(Tr./Tc).^6)) + ...            
            (exp(D + E.*(Tr./Tc)))./(pc.*vc).*(pr - pLV));            
        rhoH2 = const.MH2./vH2;                                             % density of H2 in liquid phase (g l^-1)
        
        c = const.vcH2*exp(-2.1288 - 0.1880.*(pr./pc) + ...                 % effect of composition (l mol^-1)
            5.4822*(Tr./Tc).^2);

        rho = (x18*const.M18 + x0*const.M0 + xH2*const.MH2)./ ...           % density of liquid mixture, containing H2 (g l^-1)
            ((x18 + x0).*v + xH2.*vH2 + c.*xH2.^2);
    end
    
    %% ideal gas law
    if nargout >= 3  
        ZId = 1;                                                            % compressibility factor of ideal gas mixture (-)
        y18 = x18.*gamma18.*pLV18./pr;                                      % molar fraction of H18-DBT in gas (vapour) phase (-)
        y0 = x0.*gamma0.*pLV0./pr;                                          % molar fraction of H0-DBT in gas (vapour) phase (-)
        yH2 = 1 - y18 - y0;                                                 % molar fraction of H2 in gas (vapour) phase (-)
        
        MId = y18*const.M18 + y0*const.M0 + yH2*const.MH2;                  % molar mass of ideal gas mixture (g mol^-1)
        rhoId = ((pr*1e+05).*MId./(ZId*const.R.*Tr))/1000;                  % density of ideal gas mixture (g l^-1)
        %% real gas (Soave-Redlich-Kwong EOS)
        if nargout == 4
            vle = VLE(const, pr, Tr, HG, 1);      
            ZG = vle.ZG;                                                    % compressibility factor of real gas mixture (-)    
            y18 = vle.y18;                                                  % molar fraction of H18-DBT in gas (vapour) phase (-) 
            y0 = vle.y0;                                                    % molar fraction of H0-DBT in gas (vapour) phase (-)
            yH2 = vle.yH2;                                                  % molar fraction of H2 in gas (vapour) phase (-)

            MG = y18*const.M18 + y0*const.M0 + yH2*const.MH2;               % molar mass of real gas mixture (g mol^-1)
            rhoG = ((pr*1e+05).*MG./(ZG*const.R.*Tr))/1000;                 % density of real gas mixture (g l^-1)
        end
    end
end