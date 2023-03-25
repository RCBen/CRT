function [Poy18, Poy0, PoyH2] = getPoynting(const, pr, Tr, HG)
    % [Poy18, Poy0, PoyH2] = getPoynting(const, pr, Tr, HG)
    %
    % Poynting factor of H18-DBT (-)
    %                    H0-DBT  (-)
    %                    H2      (-)
    %
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (p. 592)
    %    ISBN: 0-13-977745-8 (p. 608)
    
    nP = numel(pr);
    nTheta = numel(Tr);
    if ~isrow(pr)
        error('pr has to be a scalar or a row-vector')
    elseif isrow(pr) && (nP ~= nTheta || ~isscalar(pr) && pr(1) ~= pr(2))
        pr = repelem(pr,1,nTheta);
        Tr = repmat(Tr,1,nP);
    end
    if ~isrow(Tr)
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr))
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end
        
    [pLV18, pLV0] = getPLV(const, Tr);                                      % vapour pressure of H18-DBT & H0-DBT (bar)
    
    % see DOI: 10.1007/978-3-642-19981-3 (p. 592)
    v = [const.M18; const.M0]./getRho(const, [], Tr, [1;0]);                % molar volume of H18-DBT & H0-DBT (l mol^-1)
    Poy18 = exp(v(1,:)/1000.*(pr - pLV18)*1e+05./(const.R*Tr));             % Poynting factor of H18-DBT (-)
    Poy0 = exp(v(2,:)/1000.*(pr - pLV0)*1e+05./(const.R*Tr));               % Poynting factor of H0-DBT (-)
    
    % see ISBN: 0-13-977745-8 (p. 608)
    if nargout == 3
        [~, rhoH2] = getRho(const, pr, Tr, HG);                             % density of H2 in liquid phase (g l^-1) 
        vH2 = const.MH2./rhoH2;                                             % molar volume of H2 in liquid phase (l mol^-1)
        PoyH2 = exp(vH2/1000.*(pr - const.p0)*1e+05./(const.R*Tr));         % Poynting factor of H2 (-)
    end
end