function eqr = VLEreact(const, pr, Tr, real)
    % eqr = VLEreact(const, pr, Tr, [real])
    %
    % reaction pressure (bar)
    % reaction temperature (°C)
    % extent of reaction (mol)
    % equilibrium conversion of H18-DBT (-)
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
    % molar fraction of H18-DBT (-) in liquid phase 
    %                   H0-DBT  (-)
    %                   H2      (-)
    % total molar fraction in liquid phase (-)
    % molar fraction of H18-DBT (-) in gas (vapour) phase 
    %                   H0-DBT  (-)
    %                   H2      (-)   
    % total molar fraction in gas (vapour) phase (-) 
    % number of iteration steps needed (-)
    % solution table (-)
    % 
    % pr (bar)          scalar / vector
    % Tr (K)            scalar / vector
    
    if nargin <= 3
        real = 0;
    end
    
    nP = numel(pr);
    nTheta = numel(Tr);
    if ~isvector(pr) || ~isvector(Tr)
        error('pr & Tr have to be vectors')
    end
    
    optionsLSQ = optimoptions('lsqnonlin', ...
        'Display', 'off', 'FiniteDifferenceType', 'central', ...
        'MaxFunctionEvaluations', 1000);
    
    epsilon = []; ...                                                       % molar fraction of liquid phase (-)
        lambda = []; ...                                                    % extent of reaction (mol)
        zeq = []; ...                                                       % overall molar fraction of H18-DBT & H0-DBT & H2 (-)
        y18eq = []; y0eq = []; yH2eq = []; ...                              % molar fraction of H18-DBT & H0-DBT & H2 in gas (vapour) phase (-)
        y180eq = []; y00eq = []; yH20eq = []; ...
        x18eq = []; x0eq = []; xH2eq = []; ...                              % molar fraction of H18-DBT & H0-DBT & H2 in liquid phase (-)        
        x180eq = []; x00eq = []; xH20eq = [];
    
    n0 = [1; 0; 0];                                                         % initial molar amount of H18-DBT & H0-DBT & H2 (mol)
    [pLV18, pLV0, fPH2] = getPLV(const, Tr);                                % vapour pressure of H18-DBT & H0-DBT & Prausnitz-fugacity of H2 (bar)
    pLV = [pLV18; pLV0; fPH2];
    K = getK(const, pr, Tr);                                                % equilibrium constant for dehydrogenation of H18-DBT (-)
    
    gamma18 = ones(nP,nTheta);                                              % activity coefficient of H18-DBT (-)
    gamma0 = ones(nP,nTheta);                                               % activity coefficient of H0-DBT (-)
    gammaH2 = ones(nP,nTheta);                                              % activity coefficient of H2 (-)
    phi18V = ones(nP,nTheta);                                               % fugacity coefficient of H18-DBT @ pr (gas mixture) (-)
    phi0V = ones(nP,nTheta);                                                % fugacity coefficient of H0-DBT @ pr (gas mixture) (-)
    phiH2G = ones(nP,nTheta);                                               % fugacity coefficient of H2 @ pr (gas mixture) 
    phi18LV = ones(nP,nTheta);                                              % fugacity coefficient of H18-DBT @ pLV18 (saturation) (-)
    phi0LV = ones(nP,nTheta);                                               % fugacity coefficient of H0-DBT @ pLV0 (saturation) (-)
    phiH2LV = ones(nP,nTheta);                                              % fugacity coefficient of H2 @ fPH2 (saturation) (-)
    Poy18 = ones(nP,nTheta);                                                % Poynting factor of H18-DBT (-)
    Poy0 = ones(nP,nTheta);                                                 % Poynting factor of H0-DBT (-)
    PoyH2 = ones(nP,nTheta);                                                % Poynting factor of H2 (-)
    
    x = 0.5*ones(nP, nTheta, 2);                                            % initial point of epsilon & lambda
    lb = zeros(1,1,2);                                                      % lower bounds of epsilon & lambda
    ub = ones(1,1,2);                                                       % upper bounds of epsilon & lambda
    
    RowNames = cell(nP, 1);
    VariableNames = cell(1, nTheta);
    Nstep = zeros(nP, nTheta);                                              % number of iteration steps (-)
    residual = zeros(nP, nTheta, 3);
    exitflag = zeros(nP, nTheta);
    
    maxIter = 50;
    step = 1;
    delta = 1e-04;
    
    while step <= maxIter
        for i = 1:nP
            RowNames(i) = cellstr(strcat('p', num2str(round(pr(i), 1))));
            for j = 1:nTheta
                VariableNames(j) = cellstr(strcat('T', ...
                    num2str(round(Tr(j) - 273.15, 0))));               
                if step <= 2 || sol(i,j) == 0
                    Nstep(i,j) = step;
                    [x(i,j,:), ~, residual(i,j,:), exitflag(i,j)] = ...
                        lsqnonlin(@(x) VLEreactFun(x, const), x(i,j,:), ...
                        lb, ub, optionsLSQ);
                end
            end
        end
        
        sol = ones(nP, nTheta);                                             % auxiliary variable for defining solution condition to exit while loop
        sol(exitflag < 0 | any(zeq < -delta, 3) | ...                       %   solution was found (sol = 1) or not (sol = 0) 
            any(abs(residual) > delta, 3)) = 0;    
        
        if ~isempty(x180eq) && numel(x180eq) == numel(x18eq)
            sol(abs(x180eq - x18eq) > delta | ...
                abs(x00eq - x0eq) > delta | ...
                abs(xH20eq - xH2eq) > delta | ...
                abs(y180eq - y18eq) > delta | ...
                abs(y00eq - y0eq) > delta | ...
                abs(yH20eq - yH2eq) > delta) = 0;
        end

        if all(Nstep(:)>=2) && all(sol(:) == 1)              
            break
        else
            step = step + 1;                       
        end
    end
       
    eqr.pr = pr;
    eqr.Tr = Tr;
    eqr.epsilon = epsilon;                                                  % molar fraction of liquid phase (-)
    eqr.lambda = lambda;                                                    % extent of reaction (mol)
    eqr.X = eqr.lambda/n0(1);                                               % equilibrium conversion of H18-DBT (-)
    eqr.gamma18 = gamma18;                                                  % activity coefficient of H18-DBT (-)
    eqr.gamma0 = gamma0;                                                    % activity coefficient of H0-DBT (-)
    eqr.gammaH2 = gammaH2;                                                  % activity coefficient of H2 (-)
    eqr.phi18V = phi18V;                                                    % fugacity coefficient of H18-DBT @ pr (gas mixture) (-)
    eqr.phi0V = phi0V;                                                      % fugacity coefficient of HO-DBT @ pr (gas mixture) (-)
    eqr.phiH2G = phiH2G;                                                    % fugacity coefficient of H2 @ pr (gas mixture) (-)
    eqr.phi18LV = phi18LV;                                                  % fugacity coefficient of H18-DBT @ pLV18 (saturation) (-)
    eqr.phi0LV = phi0LV;                                                    % fugacity coefficient of H0-DBT @ pLV0 (saturation) (-)
    eqr.phiH2LV = phiH2LV;                                                  % fugacity coefficient of H2 @ fPH2 (saturation) (-)
    eqr.Poy18 = Poy18;                                                      % Poynting factor of H18-DBT (-)
    eqr.Poy0 = Poy0;                                                        % Poynting factor of H0-DBT (-)
    eqr.PoyH2 = PoyH2;                                                      % Poynting factor of H2 (-)
    eqr.x18 = x18eq;                                                        % molar fraction of H18-DBT in liquid phase @ equilibrium (-)
    eqr.x0 = x0eq;                                                          % molar fraction of H0-DBT in liquid phase @ equilibrium (-)
    eqr.xH2 = xH2eq;                                                        % molar fraction of H2 in liquid phase @ equilibrium (-)
    eqr.xges = x18eq + x0eq + xH2eq;
    eqr.y18 = y18eq;                                                        % molar fraction of H18-DBT in gas (vapour) phase @ equilibrium (-)
    eqr.y0 = y0eq;                                                          % molar fraction of H0-DBT in gas (vapour) phase @ equilibrium (-)
    eqr.yH2 = yH2eq;                                                        % molar fraction of H2 in gas (vapour) phase @ equilibrium (-)
    eqr.yges = y18eq + y0eq + yH2eq;
    eqr.Nstep = Nstep;
    eqr.tbl = array2table(sol, 'RowNames', RowNames, 'VariableNames', ...
        VariableNames);    
    
    function y = VLEreactFun(x, const)
        epsilon(i,j) = x(1);
        lambda(i,j) = x(2);

        if step > 1
            [gamma18(i,j), gamma0(i,j), gammaH2(i,j)] = ...                 
                RST(const, Tr(j), x18eq(i,j), x0eq(i,j));    
            if real == 1
                HG = x18eq(i,j)/(1-xH2eq(i,j));            
                [~, phi18V(i,j), phi0V(i,j), phiH2G(i,j), ...
                    phi18LV(i,j), phi0LV(i,j), phiH2LV(i,j)] = ...
                    SRK(const, pr(i), Tr(j), y18eq(i,j), y0eq(i,j));
                [Poy18(i,j), Poy0(i,j), PoyH2(i,j)] = ...
                    getPoynting(const, pr(i), Tr(j), HG);
            end
            
            x180eq(i,j) = x18eq(i,j); 
            x00eq(i,j) = x0eq(i,j);
            xH20eq(i,j) = xH2eq(i,j);
            y180eq(i,j) = y18eq(i,j);
            y00eq(i,j) = y0eq(i,j);
            yH20eq(i,j) = yH2eq(i,j);
        end
        
        gamma = [gamma18(i,j); gamma0(i,j); gammaH2(i,j)];                  % activity coefficient of H18-DBT & H0-DBT & H2 (-)
        phiV = [phi18V(i,j); phi0V(i,j); phiH2G(i,j)];                      % fugacity coefficients @ pr (gas mixture) (-)
        phiLV = [phi18LV(i,j); phi0LV(i,j); phiH2LV(i,j)];                  % fugacity coefficients @ pLV (saturation) (-)
        Poy = [Poy18(i,j); Poy0(i,j); PoyH2(i,j)];                          % Poynting factors (-)
        F = phiLV./phiV.*Poy;
        
        zeq(i,j, :) = (n0 + const.nu*lambda(i,j))./ ...                     % molar fraction of H18-DBT & H0-DBT & H2 in liquid & gas (vapour) phase @ equilibrium (-)
            sum(n0 + const.nu*lambda(i,j));
        
        xeq = reshape(zeq(i,j, :), 3, 1)*pr(i)./(epsilon(i,j)*pr(i) + ...
            (1-epsilon(i,j))*gamma.*pLV(:,j).*F);        
        x18eq(i,j) = xeq(1);
        x0eq(i,j) = xeq(2);
        xH2eq(i,j) = xeq(3);
        
        yeq = (reshape(zeq(i,j, :), 3, 1) - xeq*epsilon(i,j))/ ...
            (1-epsilon(i,j));
        y18eq(i,j) = yeq(1);
        y0eq(i,j) = yeq(2);
        yH2eq(i,j) = yeq(3);        
        
        preq = sum(xeq.*gamma.*pLV(:,j).*F);                                % reaction pressure @ equilibrium (bar)
        
        a18 = x18eq(i,j)*gamma(1);                                          % activity of H18-DBT (-)
        a0 =  x0eq(i,j)*gamma(2);                                           % activity of H0-DBT (-)
        fH2 = yH2eq(i,j)*pr(i)/const.p0*phiV(3);                            % fugacity of H2 (-)
        
        y = [preq/pr(i) - 1; ...
            log(prod([a18; a0; fH2].^const.nu)/K(j)); ...                   % K = a0*fH2^9/a18 = (gamma0*c0)*(phiH2G*pH2/p0)^9/(gamma18*c18) (-)                 
            x18eq(i,j) + x0eq(i,j) + xH2eq(i,j) - 1];                       % additional equation in order to prevent any(x) > 1
    end
end