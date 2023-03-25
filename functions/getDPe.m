function [DPe18, DPe0, DPeH2] = getDPe(const, p, Tr, HG)
    % [DPe18, DPe0, DPeH2] = getDPe(const, p, Tr, HG)
    %
    % binary diffusion coefficient of H18-DBT (m^2 s^-1) in liquid mixture, wo H2 in catalyst pellet
    %                                 H0-DBT  (m^2 s^-1) 
    %                                 H2      (m^2 s^-1)
    % 
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (p. 171 && p.653)
    %          10.1002/aic.690010222
    %    ISBN: 0-07-0011682-2 (pp. 11.21 f & 11.23 ff & 11.35 & 11.42)
    
    nTheta = numel(Tr);
    if ~isrow(Tr) 
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end
    
    F = 1;                                                                  % restrictive diffusion factor (-)
    if p.FPe == 1
        if ~isfield(p, 'dpore')
            p.dpore = const.dpore;                                          % average pore diameter (nm) 
        end  
        F = (1-const.dmol18./p.dpore).^4.9;                                 
    end
    
    %% Tyn and Calus method (x < 0.05 - 0.10)
    % see DOI: 10.1007/978-3-642-19981-3 (p. 171)
    %    ISBN: 0-07-0011682-2 (pp. 11.23 ff)
    [~, P18, P0] = getSigma(const, Tr, HG);                                 % Parachor of H18-DBT & H0-DBT (cm^3 g^0.25 s^-0.5 mol^-1)
    eta  = getEta(const, [], Tr, [1;0]);                                    % dynamic viscosity of H18-DBT & H0-DBT (mPas)    
    D180inf = (8.93e-08*(1000*const.vB18/(1000*const.vB0)^2)^(1/6)* ...     % binary diffusion coefficient of H18-DBT in H0-DBT @ infinite dilution (m^2 s^-1)
        (P0./P18).^0.6.*Tr./eta(2,:))*1e-04;                                 
    D018inf = (8.93e-08*(1000*const.vB0/(1000*const.vB18)^2)^(1/6)* ...     % binary diffusion coefficient of H0-DBT in H18-DBT @ infinite dilution (m^2 s^-1)
        (P18./P0).^0.6.*Tr./eta(1,:))*1e-04;                                

    %% Vignes correlation
    % see DOI: 10.1007/978-3-642-19981-3 (p. 171)
    %    ISBN: 0-07-0011682-2 (p. 11.35)
    deltaHG = 1e-06;                                                        % length of interval for calculating the 1st derivative (-)
    HGm = HG - deltaHG/2;                                                   % lower boundary (-)
    HGp = HG + deltaHG/2;                                                   % upper boundary (-)
    HGm(HGm < 0) = 0;
    HGp(HGp > 1) = 1;
    gamma18 = RST(const, Tr, [HGm; HGp]);                                   % activity coefficient of H18-DBT (-)
    alpha = (log(HGm.*gamma18(1:size(HG,1),:)) - ...                        % thermodynamic correction due to non-ideality of liquid phase (-)
        log(HGp.*gamma18(size(HG,1)+1:end,:)))./(log(HGm) - log(HGp));
    alpha(isnan(alpha)) = 1;
    D180 = (D180inf.^(1-HG).*D018inf.^HG).*alpha;                           % binary diffusion coefficient of H18-DBT in H0-DBT (m^2 s^-1)
    D018 = (D018inf.^HG.*D180inf.^(1-HG)).*alpha;                           % binary diffusion coefficient of H0-DBT in H18-DBT (m^2 s^-1)
    
    switch p.DPe
        case 1
            %1 - ratio 1/10
            fPe = 0.1;
        case 2
            %2 - Dumanski
            fPe = -const.epsPe.*((1 - const.epsPe).^(2/3) - 1);        
        case 3
            %3 - Hugo
            fPe = const.epsPe.*(5.*const.epsPe-1)./4;            
        case 4
            %4 - Bruggemann
            fPe = const.epsPe.^(3/2);
        case 5
            %5 - Maxwell
            fPe = -const.epsPe./(const.epsPe./2 - 3/2);
        case 6
            %6 - Millington/Quirk
            fPe = const.epsPe.^(4/3);
    end
    
    DPe18 = fPe.*F.*D180;                                                   % effective diffusion coefficient of H18-DBT in liquid mixture, wo H2 (m^2 s^-1)
    DPe0 = fPe.*F.*D018;                                                    % effective diffusion coefficient of H0-DBT in liquid mixture, wo H2 (m^2 s^-1)
    
    %% Wilke-Chang estimation method (x < 0.05 - 0.10)
    % see DOI: 10.1002/aic.690010222
    %    ISBN: 0-07-0011682-2 (pp. 11.21 f & 11.42)
    if nargout == 3
        phi18 = 1;                                                          % association factor of H18-DBT (-)
        phi0 = 1;                                                           % association factor of H0-DBT (-)
        DH2inf = (7.4e-08.* ...                                             % binary diffusion coefficient of H2 in liquid mixture, wo H2 @ infinite dilution (m^2 s^-1)
            sqrt(HG*phi18*const.M18 + (1-HG)*phi0*const.M0).*Tr./ ... 
            (getEta(const, [], Tr, HG)*14.3^0.6))*1e-04;   

        DPeH2 = fPe.*F.*DH2inf;                                             % effective diffusion coefficient of H2 in liquid mixture, wo H2 (m^2 s^-1)
    end
end