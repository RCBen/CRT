function [lambda, lambdaH2Id, lambdaId, lambdaG] = getLambda(const, pr, Tr, HG)
    % [lambda, lambdaH2Id, lambdaId, lambdaG] = getLambda(const, pr, Tr, HG)
    %
    % heat conductivity of liquid mixture, wo H2 (W m^-1 K^-1)
    %                      ideal gas H2          (W m^-1 K^-1)
    %                      ideal gas mixture     (W m^-1 K^-1)
    %                      real gas mixture      (W m^-1 K^-1)
    %
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 166 ff & 357 & 445 & 554 f)
    %    ISBN: 0-07-0011682-2 (pp. 10.12 & 10.18 & 10.30 f & 10.45 f & 10.57)
    
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
    if ~isrow(Tr)
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end
    
    %% Sastri method
    % see ISBN: 0-07-0011682-2 (pp. 10.45 f) 
    lambdaB18 = 1*0.0545 + (2+13)*-0.0008 + 5*-0.0600 + 3*0.1130;
    lambda18 = lambdaB18* ...                                               % heat conductivity of H18-DBT (W m^-1 K^-1)
        0.16.^(1 - ((1 - Tr/const.Tc18)/(1 - const.Tb18/const.Tc18)).^0.2);
  
    %% linear fit to VDI data
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 554 f)
    lambda0 = 0.1693 - 0.0001318*Tr;                                        % heat conductivity of H0-DBT (W m^-1 K^-1)
        
    %% Jamieson correlation
    % see ISBN: 0-07-0011682-2 (p. 10.57) 
    w18 = HG*const.M18./(HG*const.M18 + (1-HG)*const.M0);                   % mass fraction of H18-DBT in liquid phase, wo H2 (-)
    w0 = 1 - w18;                                                           % mass fraction of H0-DBT in liquid phase, wo H2 (-)
%     lambda = w18.*lambda18 + w0.*lambda0 - ...                              % heat conductivity of liquid mixture, wo H2 (W m^-1 K^-1)
%         0.72*w18.*w0.*abs(lambda0 - lambda18);    
    lambda = w18.*lambda18 + w0.*lambda0 - ...                              % heat conductivity of liquid mixture, wo H2 (W m^-1 K^-1)
        1*w18.*w0.*abs(lambda0 - lambda18).*(1-sqrt(w0)).*w0;
    
    %% VDI data
    % see DOI: 10.1007/978-3-642-19981-3 (p. 357 (eq. 5) & 445)
    lambdaH2Id = (0.651e-03 + 0.76730e-03*Tr - 0.687050e-06*Tr.^2 + ...     % ideal gas heat conductivity of H2 (W m^-1 K^-1)
        0.506510e-09*Tr.^3 - 0.138540e-12*Tr.^4);  
       
    %% method of Chung
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 167 f)
    %    ISBN: 0-07-0011682-2 (p. 10.12)
    if nargout >= 3
        [~, ~, etaId]  = getEta(const, pr, Tr, [1;0]);                      % ideal gas dynamic viscosity of H18-DBT & H0-DBT (uPas)
        [~, ~, cpId] = getCp(const, pr, Tr, [1;0]);                         % ideal gas heat capacity of H18-DBT & H0-DBT @ constant pressure (kJ mol^-1 K^-1)
        
        alpha18 = cpId(1,:)*1000/const.R - 2.5;
        beta18 = 0.7862 - 0.7109*const.w18 + 1.3168*const.w18^2;
        gamma18 = 2 + 10.5*(Tr/const.Tc18).^2;        
        psi18 = 1 + alpha18.*((0.215 + 0.28288*alpha18 - 1.061*beta18 + ...
            0.26665*gamma18)./ ...
            (0.6366 + beta18*gamma18 + 1.061*alpha18*beta18));       
        lambda18Id = 3.75*psi18.*(etaId(1,:)*1e-06)*const.R/ ...            % ideal gas heat conductivity of H18-DBT (W m^-1 K^-1)
            (const.M18/1000);
             
        alpha0 = cpId(2,:)*1000/const.R - 2.5;
        beta0 = 0.7862 - 0.7109*const.w0 + 1.3168*const.w0^2;
        gamma0 = 2 + 10.5*(Tr/const.Tc0).^2;        
        psi0 = 1 + alpha0.*((0.215 + 0.28288*alpha0 - 1.061*beta0 + ...
            0.26665*gamma0)./ ...
            (0.6366 + beta0*gamma0 + 1.061*alpha0*beta0));       
        lambda0Id = 3.75*psi0.*(etaId(2,:)*1e-06)*const.R/ ...              % ideal gas heat conductivity of H0-DBT (W m^-1 K^-1)
            (const.M0/1000); 
        
        %% method of Wassiljewa, Mason & Saxena
        % see ISBN: 0-07-0011682-2 (pp. 10.30 f) 
        if size(HG, 1) == 2 && all(HG == [1;0])
            lambdaId = [lambda18Id; lambda0Id];                             % ideal gas heat conductivity of H18-DBT & H0-DBT (W m^-1 K^-1)
        else
            vle = VLE(const, pr, Tr, HG);           
            y18 = vle.y18;                                                  % molar fraction of H18-DBT in gas (vapour) phase (-) 
            y0 = vle.y0;                                                    % molar fraction of H0-DBT in gas (vapour) phase (-)
            yH2 = vle.yH2;                                                  % molar fraction of H2 in gas (vapour) phase (-)
            
            Tr18 = Tr/const.Tc18;                                           % reduced temperature of H18-DBT (K)
            Tr0 = Tr/const.Tc0;                                             % reduced temperature of H0-DBT (K)
            TrH2 = Tr/const.TcH2q;                                          % reduced temperature of H2 (K)

            Gamma18 = 210*(const.Tc18*const.M18^3*const.pc18^-4)^(1/6);     % reduced inverse thermal conductivity of H18-DBT (W m^-1 K^-1)^-1
            Gamma0 = 210*(const.Tc0*const.M0^3*const.pc0^-4)^(1/6);         % reduced inverse thermal conductivity of H0-DBT (W m^-1 K^-1)^-1
            GammaH2 = 210*(const.TcH2q*const.MH2^3*const.pcH2q^-4)^(1/6);   % reduced inverse thermal conductivity of H2 (W m^-1 K^-1)^-1

            lambda180 = ...
                (Gamma0*(exp(0.0464*Tr18) - exp(-0.2412*Tr18)))./ ...
                (Gamma18*(exp(0.0464*Tr0) - exp(-0.2412*Tr0)));
            lambda18H2 = ...
                (GammaH2*(exp(0.0464*Tr18) - exp(-0.2412*Tr18)))./ ...
                (Gamma18*(exp(0.0464*TrH2) - exp(-0.2412*TrH2)));
            lambda018 =lambda180.^-1;
            lambda0H2 = ...
                (GammaH2*(exp(0.0464*Tr0) - exp(-0.2412*Tr0)))./ ...
                (Gamma0*(exp(0.0464*TrH2) - exp(-0.2412*TrH2)));
            lambdaH218 = lambda18H2.^-1;
            lambdaH20 = lambda0H2.^-1;

            A1818 = 1;
            A180 = (1 + sqrt(lambda180)*(const.M18/const.M0)^(1/4)).^2/ ...
                sqrt(8*(1 + const.M18/const.M0));
            A18H2 = (1 + sqrt(lambda18H2)*(const.M18/const.MH2)^(1/4)).^2/ ...
                sqrt(8*(1 + const.M18/const.MH2));
            A018 = (1 + sqrt(lambda018)*(const.M0/const.M18)^(1/4)).^2/ ...
                sqrt(8*(1 + const.M0/const.M18));
            A00 = 1;
            A0H2 = (1 + sqrt(lambda0H2)*(const.M0/const.MH2)^(1/4)).^2/ ...
                sqrt(8*(1 + const.M0/const.MH2));
            AH218 = (1 + sqrt(lambdaH218)*(const.MH2/const.M18)^(1/4)).^2/ ...
                sqrt(8*(1 + const.MH2/const.M18));
            AH20 = (1 + sqrt(lambdaH20)*(const.MH2/const.M0)^(1/4)).^2/ ...
                sqrt(8*(1 + const.MH2/const.M0));
            AH2H2 = 1;

            lambdaId = ...                                                  % heat conductivity of ideal gas mixture (W m^-1 K^-1) 
                y18.*lambda18Id./(y18.*A1818 + y0.*A180 + yH2.*A18H2) + ...
                y0.*lambda0Id./(y18.*A018 + y0.*A00 + yH2.*A0H2) + ...
                yH2.*lambdaH2Id./(y18.*AH218 + y0.*AH20 + yH2.*AH2H2);       
        end
        %% neglectable effect of pressure (pr < 10 bar) 
        % see ISBN: 0-07-0011682-2 (p. 10.18)  
        if nargout == 4
            lambdaG = lambdaId;                                             % heat conductivity of real gas mixture (W m^-1 K^-1)
        end
    end
end