function [sigma, P18, P0] = getSigma(const, Tr, HG)
    % [sigma, P18, P0] = getSigma(const, Tr, HG)
    %
    % surface tension of liquid mixture, wo H2 (mN m^-1)
    % Parachor of H18-DBT (cm^3 g^0.25 s^-0.5 mol^-1)
    %             H0-DBT  (cm^3 g^0.25 s^-0.5 mol^-1)
    %
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 169 f)
    %    ISBN: 0-07-0011682-2 (pp. 12.2 ff & 12.13 f) 
    
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

    %% Brock/Bird/Miller corresponding states correlation
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 169 f)
    %    ISBN: 0-07-0011682-2 (pp. 12.3 ff)
    Q18 = 0.1196*(1 + (const.Tb18/const.Tc18*log(const.pc18/1.01325))/ ...
        (1 - const.Tb18/const.Tc18)) - 0.279;
    sigma18 = const.pc18^(2/3)*const.Tc18^(1/3)*Q18* ...                    % surface tension of H18-DBT (mN m^-1)
        (1 - Tr/const.Tc18).^(11/9);
    
    Q0 = 0.1196*(1 + (const.Tb0/const.Tc0*log(const.pc0/1.01325))/ ...
        (1 - const.Tb0/const.Tc0)) - 0.279;
    sigma0 = const.pc0^(2/3)*const.Tc0^(1/3)*Q0* ...                        % surface tension of H0-DBT (mN m^-1)
        (1 - Tr/const.Tc0).^(11/9);
 
    %% Macleod-Sugden correlation
    % see DOI: 10.1007/978-3-642-19981-3 (p. 170)
    %    ISBN: 0-07-0011682-2 (pp. 12.13 f)
    v = [const.M18; const.M0]./(getRho(const, [], Tr, [1;0])/1000);         % molar volume of H18-DBT & H0-DBT (cm^3 mol^-1)
    P18 = v(1,:).*sigma18.^(1/4);                                           % Parachor of H18-DBT (cm^3 g^0.25 s^-0.5 mol^-1)
    P0 = v(2,:).*sigma0.^(1/4);                                             % Parachor of H0-DBT (cm^3 g^0.25 s^-0.5 mol^-1)
    P = HG.^2.*P18 + 2*HG.*(1-HG).*(P18 + P0)/2 + (1-HG).^2.*P0;            % Parachor of liquid mixture, wo H2 (cm^3 g^0.25 s^-0.5 mol^-1)
    sigma = (P./((HG*const.M18 + (1-HG)*const.M0)./ ...                     % surface tension of liquid mixture, wo H2 (mN m^-1)
        (getRho(const, [], Tr, HG)/1000))).^3.6;                                                     
end