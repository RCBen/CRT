function [DEff18, DEff0, DEffH2] = getDEff(const, p, Tr, HG)
    % [DEff18, DEff0, DEffH2] = getDEff(const, p, Tr, HG)
    %
    % binary diffusion coefficient of H18-DBT (m^2 s^-1) in liquid mixture, wo H2 in porous domain (bed)
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
    
    switch p.DEff
        case 1
            %1 - ratio 1/10
            fEff = 0.1;
        case 2
            %2 - Dumanski
            fEff = -p.epsBed.*((1 - p.epsBed).^(2/3) - 1);
        case 3
            %3 - Hugo
            fEff = p.epsBed.*(5.*p.epsBed-1)./4;          
        case 4
            %4 - Bruggemann
            fEff = p.epsBed.^(3/2);
        case 5
            %5 - Maxwell
            fEff = -p.epsBed./(p.epsBed./2 - 3/2);
        case 6 
            %6 - Millington/Quirk
            fEff = p.epsBed.^(4/3);
    end
    
    DEff18 = fEff.*D180;                                                    % effective diffusion coefficient of H18-DBT in liquid mixture, wo H2 (m^2 s^-1)
    DEff0 = fEff.*D018;                                                     % effective diffusion coefficient of H0-DBT in liquid mixture, wo H2 (m^2 s^-1)
    
    %% Wilke-Chang estimation method (x < 0.05 - 0.10)
    % see DOI: 10.1002/aic.690010222
    %    ISBN: 0-07-0011682-2 (pp. 11.21 f & 11.42)
    if nargout == 3
        phi18 = 1;                                                          % association factor of H18-DBT (-)
        phi0 = 1;                                                           % association factor of H0-DBT (-)
        DH2inf = (7.4e-08.* ...                                             % binary diffusion coefficient of H2 in liquid mixture, wo H2 @ infinite dilution (m^2 s^-1)
            sqrt(HG*phi18*const.M18 + (1-HG)*phi0*const.M0).*Tr./ ... 
            (getEta(const, [], Tr, HG)*14.3^0.6))*1e-04;   

        DEffH2 = fEff*DH2inf;                                               % effective diffusion coefficient of H2 in liquid mixture, wo H2 (m^2 s^-1)
    end
end