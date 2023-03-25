function vle = VLE(const, pr, Tr, HG, real)
    % vle = VLE(const, pr, Tr, HG, [real])
    %
    % compressibility factor of gas mixture (-)
    % molar fraction of H18-DBT (-) in gas (vapour) phase
    %                   H0-DBT  (-)
    %                   H2      (-)
    % molar fraction of H18-DBT (-) in liquid phase
    %                   H0-DBT  (-)
    %                   H2      (-)
    % activity coefficient of H18-DBT (-)
    %                         H0-DBT  (-)
    %                         H2      (-)
    % fugacity coefficient of H18-DBT (-) @ pr (V/G)
    %                         H0-DBT  (-)
    %                         H2      (-)
    % fugacity coefficient of H18-DBT (-) @ saturation (LV)
    %                         H0-DBT  (-)
    %                         H2      (-)
    % Poynting factor of H18-DBT (-)
    %                    H0-DBT  (-)
    %                    H2      (-)
    % reason the solver stopped (-)
    %
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    
    if nargin <= 4
        real = 0;
    end
    
    nP = numel(pr);
    nTheta = numel(Tr);
    if ~isrow(pr)
        error('pr has to be a scalar or a row-vector')
    elseif isscalar(pr)
        pr = repelem(pr,1,nTheta);
    elseif isrow(pr) && (nP ~= nTheta || pr(1) ~= pr(2))
        pr = repelem(pr,1,nTheta);
        Tr = repmat(Tr,1,nP);
        nTheta = numel(Tr);
    end
    if ~isrow(Tr)
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr))
        if ~isscalar(HG) && iscolumn(HG)
            HG = repmat(HG,1,nTheta);
         elseif ~isscalar(HG) && isrow(HG)
            HG = repmat(HG',1,nTheta);
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end
    
    optionsLSQ = optimoptions('lsqnonlin', ...
        'Display', 'off', 'FiniteDifferenceType', 'central', ...
        'MaxFunctionEvaluations', 1000);
    
    ZG = ones(size(HG,1),nTheta);                                           % compressibility factor of gas mixture (-)
    y18 = []; y0 = []; yH2 = []; ...                                        % molar fraction of H18-DBT & H0-DBT & H2 in gas (vapour) phase (-)
        x18 = []; x0 = []; xH2 = [];                                        % molar fraction of H18-DBT & H0-DBT & H2 in liquid phase (-)
    
    [pLV18, pLV0, fPH2] = getPLV(const, Tr);                                % vapour pressure of H18-DBT & H0-DBT & Prausnitz-fugacity of H2 (bar)
    pLV = [pLV18; pLV0; fPH2];
    [gamma18, gamma0, gammaH2] = RST(const, Tr, HG);                        % activity coefficient of H18-DBT & H0-DBT & H2 (-)
    phi18V = ones(size(HG,1),nTheta);                                       % fugacity coefficient of H18-DBT @ pr (gas mixture) (-)
    phi0V = ones(size(HG,1),nTheta);                                        % fugacity coefficient of H0-DBT @ pr (gas mixture) (-)
    phiH2G = ones(size(HG,1),nTheta);                                       % fugacity coefficient of H2 @ pr (gas mixture) 
    phi18LV = ones(size(HG,1),nTheta);                                      % fugacity coefficient of H18-DBT @ pLV18 (saturation) (-)
    phi0LV = ones(size(HG,1),nTheta);                                       % fugacity coefficient of H0-DBT @ pLV0 (saturation) (-)
    phiH2LV = ones(size(HG,1),nTheta);                                      % fugacity coefficient of H2 @ fPH2 (saturation) (-)
    switch real
        case 0
            Poy18 = ones(1,nTheta);                                         % Poynting factor of H18-DBT (-)
            Poy0 = ones(1,nTheta);                                          % Poynting factor of H0-DBT (-)
            PoyH2 = ones(size(HG,1),nTheta);                                % Poynting factor of H2 (-)
        case 1
            [Poy18, Poy0, PoyH2] = getPoynting(const, pr, Tr, HG);          % Poynting factor of H18-DBT & H0-DBT & H2 (-)
    end
    
    X(:,:,1) = HG.*gamma18.*pLV18./pr;                                      % initial point of y18
    X(:,:,2) = (1-HG).*gamma0.*pLV0./pr;                                    % initial point of y0
    lb = zeros(1,1,2);                                                      % lower bounds of y18& y0
    ub = ones(1,1,2);                                                       % upper bounds of y18& y0
    
    residual = zeros(size(X));
    exitflag = zeros(size(HG,1), nTheta);
    
    for i = 1:size(HG,1)
        for j = 1:nTheta
            [X(i,j,:), ~, residual(i,j,:), exitflag(i,j)] = ...
                lsqnonlin(@(x) VLEFun(x, const), X(i,j,:), lb, ub, ...
                optionsLSQ);
        end
    end
    
    vle.ZG = ZG;                                                            % compressibility factor of gas mixture (-)
    vle.y18 = y18;                                                          % molar fraction of H18-DBT in gas (vapour) phase (-)
    vle.y0 = y0;                                                            % molar fraction of H0-DBT in gas (vapour) phase (-)
    vle.yH2 = yH2;                                                          % molar fraction of H2 in gas (vapour) phase (-)
    vle.x18 = x18;                                                          % molar fraction of H18-DBT in liquid phase (-)
    vle.x0 = x0;                                                            % molar fraction of H0-DBT in liquid phase (-)
    vle.xH2 = xH2;                                                          % molar fraction of H2 in liquid phase (-)
    vle.gamma18 = gamma18;                                                  % activity coefficient of H18-DBT (-)
    vle.gamma0 = gamma0;                                                    % activity coefficient of H0-DBT (-)
    vle.gammaH2 = gammaH2;                                                  % activity coefficient of H2 (-)
    vle.phi18V = phi18V;                                                    % fugacity coefficient of H18-DBT @ pr (gas mixture) (-)
    vle.phi0V = phi0V;                                                      % fugacity coefficient of HO-DBT @ pr (gas mixture) (-)
    vle.phiH2G = phiH2G;                                                    % fugacity coefficient of H2 @ pr (gas mixture) (-)
    vle.phi18LV = phi18LV;                                                  % fugacity coefficient of H18-DBT @ pLV18 (saturation) (-)
    vle.phi0LV = phi0LV;                                                    % fugacity coefficient of H0-DBT @ pLV0 (saturation) (-)
    vle.phiH2LV = phiH2LV;                                                  % fugacity coefficient of H2 @ fPH2 (saturation) (-)
    vle.Poy18 = Poy18;                                                      % Poynting factor of H18-DBT (-)
    vle.Poy0 = Poy0;                                                        % Poynting factor of H0-DBT (-)
    vle.PoyH2 = PoyH2;                                                      % Poynting factor of H2 (-)
    vle.exitflag = exitflag;
    
    function Y = VLEFun(X, const)
        y18(i,j) = X(1);                                                    % molar fraction of H18-DBT in gas (vapour) phase (-)
        y0(i,j) = X(2);                                                     % molar fraction of H0-DBT in gas (vapour) phase (-)
        yH2(i,j) = 1 - y18(i,j) - y0(i,j);                                  % molar fraction of H2 in gas (vapour) phase (-)
        
        if numel(x18) == (i-1)*nTheta+j-1
            x18(i,j) = HG(i,j);
            x0(i,j) = 1-HG(i,j);
        else
            [gamma18(i,j), gamma0(i,j), gammaH2(i,j)] = ...                 
                RST(const, Tr(j), x18(i,j), x0(i,j)); 
        end
        if real == 1
            [ZG(i,j), phi18V(i,j), phi0V(i,j), phiH2G(i,j), ...
                phi18LV(i,j), phi0LV(i,j), phiH2LV(i,j)] = ...
                SRK(const, pr(j), Tr(j), y18(i,j), y0(i,j));
        end
        
        gamma = [gamma18(i,j); gamma0(i,j); gammaH2(i,j)];                  % activity coefficient of H18-DBT & H0-DBT & H2 (-)
        phiV = [phi18V(i,j); phi0V(i,j); phiH2G(i,j)];                      % fugacity coefficients @ pr (gas mixture) (-)
        phiLV = [phi18LV(i,j); phi0LV(i,j); phiH2LV(i,j)];                  % fugacity coefficients @ pLV (saturation) (-)
        Poy = [Poy18(1,j); Poy0(1,j); PoyH2(i,j)];                          % Poynting factors (-)
        F = phiLV./phiV.*Poy;
        
        xH2(i,j) = (pr(j) + ...                                             % molar fraction of H2 in liquid phase (-)
            sum([-HG(i,j); HG(i,j)-1].*gamma(1:2).*pLV(1:2,j).*F(1:2)))./ ...
            sum([-HG(i,j); HG(i,j)-1; 1].*gamma.*pLV(:,j).*F);
        x18(i,j) = HG(i,j)*(1-xH2(i,j));                                    % molar fraction of H18-DBT in liquid phase (-)
        x0(i,j) = (1-HG(i,j))*(1-xH2(i,j));                                 % molar fraction of H0-DBT in liquid phase (-)
        x = [x18(i,j); x0(i,j); xH2(i,j)];
        
        y = x.*gamma.*pLV(:,j).*F/pr(j);
        y18(i,j) = y(1);                                                    % molar fraction of H18-DBT in gas (vapour) phase (-)
        y0(i,j) = y(2);                                                     % molar fraction of H0-DBT in gas (vapour) phase (-)
        yH2(i,j) = y(3);                                                    % molar fraction of H2 in gas (vapour) phase (-)                          

        Y = [y18(i,j); y0(i,j)] - reshape(X,2,1);
    end
end