function [mdl, d, s, r, mb, p] = myFitLSQnl(rate, ir, irT, film, ps)
    % [mdl, d, s, r, mb, p] = myFitLSQnl(rate, ir, irT, film, ps)
    % 
    % fitting the experimental data
    % 
    % rate (dependency)         0 partial density (rho) 
    %                           1 molar concentration (c)
    %                           2 weight fraction (w)
    % i(nternal) r(esistance)   0 off
    %                           1 on
    % irT (ir regarding T)      0 off
    %                           1 on                                                                        
    % film (film limitation)    0 off
    %                           1 on
    % p(attern)s(earch)         0 off
    %                           1 on
    
    if ~exist('p', 'var') || ~isfield(p,'subdir')
        mydir = pwd;
        idcs = strfind(mydir,'\');
        p.subdir = mydir(1:idcs(end-2)-1);
        addpath(strcat(p.subdir,'\functions'))
    end
    
    %% Definition                                                             
    p.plotPS = 1;                                                           % plot patternsearch output     0 off
                                                                            %                               1 on 
    p.plot = 0;                                                             % plot results                  0 off
                                                                            %                               1 on
    p = myDefinitions(p);                                                   % includes p.rate, p.ir, p.irT, p.film, p.iso, p.DPe, p.FPe, p.DEff, p.lambdaPe, p.lambdaEff                                                                     
    if nargin <= 1
      p.ps = 0;                                                             % patternsearch                 0 off     
                                                                            %                               1 on 
    else
      p.rate = rate;
      p.ir = ir;
      p.irT = irT;
      p.film = film;
      p.ps = ps;
    end                                                                         
    
    % constants
    const = getConstants;
    
    % variables
    p.tr = 3600;                                                            % considered reaction time (s)
    p.pr = [1, 2, 3, 4, 5];                                                 % reaction pressure (bar)
    p.nP = numel(p.pr);                                                     % number of different absolute pressures (-)
    p.Tr = [310, 300, 290, 280, 270] + 273.15;                              % reaction temperature (K)
    p.nTheta = numel(p.Tr);                                                 % number of different temperatures (-)
    p.dataname = 'Data_MCK';                                                % name of file which contains experimental data
    
    % parameters
    p.Vr = [100, 600, 600, 600, 600];                                       % reactor volume (ml)
    p.Vr = p.Vr(p.pr);
    p.nL0 = [0.1, 0.65, 0.65, 0.65, 0.65];                                  % initial molar amount of liquid phase (mol);
    p.nL0 = p.nL0(p.pr);
    p.xPt = [0.0005, 0.001, 0.001, 0.001, 0.001];                           % molar fraction of Pt regarding to nL0 (-)
    p.xPt = p.xPt(p.pr);
    p.dpore = const.dpore;                                                  % average pore diameter (nm)
    
    %% Data import   
    data = cell(1,p.nP);    
    for j = 1:p.nP
        data{j} = readmatrix(p.dataname, 'Sheet', p.pr(j));
        if j == 1            
            p.I = find(data{j}(:,1)<=p.tr, 1, 'last');                      % number of considered reaction times (-)
            p.J = p.nP*p.nTheta;                                            % number of considered reactions
            p.t = data{j}(1:p.I,1);                                         % reaction time from 0 to p.tr (s)
            p.HGraw = zeros(p.I,p.J);                                       % hydrogenation grade (-)
        end
        p.HGraw(:,p.nTheta*(j-1)+1:j*p.nTheta) = data{j}(1:p.I,2:p.nTheta+1);
    end
    
    p.HG = p.HGraw;
%     if any(isnan(p.HG(:)))
%         p.HG = fillmissing(p.HGraw, 'pchip', 1, 'SamplePoints', p.t, ...
%             'EndValues', 'extrap');
%     end
%     
%     if ~isempty(find(isoutlier(p.HGraw), 1))
%         p.HG = filloutliers(p.HG, 'pchip', 1, 'SamplePoints', p.t);
%     end
    
    span = 5;                                                               % span must be odd
    for j = 1:p.J
        p.HG(:,j) = smooth(p.t, p.HG(:,j), span, 'moving');
    end
    
    %% solve mass balance
    [mb, p] = myMassBalance(p);
    p.Tm = repmat(mean(p.Tr(:)), 1, p.J);                                   % temperature in the middle of the experimental range (K)
     
%     p.urelHG = 0.4076.*p.HG.^2 - 0.7418.*p.HG + 0.3502;                     % relative uncertainty for hydrogenation grade measurements without path consideration for k = 2 (-)   
   p.urelHG = 2*(0.2641.*p.HG.^2 - 0.6199.*p.HG + 0.3884);                     % relative uncertainty for hydrogenation grade measurements considering worst case form wrong path(-) for k = 2 (-)  
%     p.urelHG = 2+(0.2038.*p.HG.^2 - 0.3709.*p.HG + 0.1751);                 % relative uncertainty for hydrogenation grade measurements without path consideration (-) for k = 2 (-) 
    p.uHG = p.urelHG.*p.HG./2;                                              % standard uncertainty for hydrogenation grade measurements (-)
    p.un18 = mb.n0r./(1 - p.HG).^2.*p.uHG;                                  % standard uncertainty for calculation of n18 (Gaussian error propagation) (-)
    I = find(p.t >= 300, 1, 'first');
    if I >= 3
        p.weight = (1./p.un18.^2).*...                                      % weighted least square fit (with response variance un18.^2) 
            repmat([log(1:(exp(1)-1)/(I-1):exp(1)), ...
            ones(1, p.I - I)]', 1, p.J);        
    end
    p.weight(1,:) = 1;                                                      % arbitrary non-zero value, because s.Residuals(1,:) = 0
    
    %% Calculation (nonlinear regression)
    Arrhenius = 'k*exp(-A*(Tm/Tr - 1))';
    dep = {'rho18', 'c18', 'w18'};
    switch p.nP
        case 1
            order = 'n1';
        otherwise 
            order = '(n1 + n2*pr)';
    end
    Rate = strcat(Arrhenius, '*', dep{p.rate+1}, '^', order);
    CoefficientNames = {'A' 'k' 'n1' 'n2'};  

    p_start = [25 5 30 -50];
    p.D = diag([1 1e-05 1e-01 1e-02]);
%     p.D = diag([1 exp(-p_start(1)) 1e-01 1e-02]);                   % equation to be solved k*exp(A.*(1-T/Tm)) with EA = A*R*Tm and k0 = k*exp(A). thus scale with exp(-A)
    if p.rate == 0
        p.D(2,2) = p.D(2,2)*1e-08;
    end

    p_lbPS = [];
    p_ubPS = [];

    switch p.ir
        case 0
            p_lb = [];
            p_ub = [];          
        case 1
%             p.D(2,2) = p.D(2,2)*1e-05;
            if p.rate == 2
                p_start = [25 20 5 100];
                p_lb = [0 0 0 0];
                p_ub = [50 1e06 1e02 5e02];
            else
%                 p_start = [40 20 30 -50];
                p_lb = [0 0 0 -5e02];
                p_ub = [50 1e06 1e02 0];
            end
%             p_lb = p_lb*p.D;
%             p_ub = p_ub*p.D;
    end

    if p.nP == 1
        CoefficientNames = CoefficientNames(1:end-1);
        p.D = p.D(1:end-1,1:end-1);
        p_start = p_start(1:end-1);
    end
    
    startTime = tic;
    
    n18rvcalc = [];
    n0rvcalc = [];
    optionsLSQ = optimoptions('lsqnonlin', 'Algorithm', ...
        'trust-region-reflective', 'Display', 'iter', ...
        'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', ...
        500*numel(p_start));
    [p_est, s.SSE, s.Residuals, exitflagLSQ, outputLSQ, ~, Jacobian] = ...
        lsqnonlin(@(p_est) myModel(p_est*p.D), p_start, ...
        p_lb, p_ub, optionsLSQ);
    calcTimeLSQ = toc(startTime);
        
    switch p.ps
        case 1
            optionsPS = optimoptions('patternsearch', 'Display', 'iter', ...   
                'AccelerateMesh', false, 'Cache', 'on', 'UseCompletePoll', true, ...              
                'SearchFcn', @MADSPositiveBasis2N, 'UseCompleteSearch', true, ...
                'UseParallel', false);            
            if p.plotPS == 1 
                optionsPS = optimoptions(optionsPS, 'PlotFcn', @psplotbestf);
            end
            [p_est, ~, exitflagPS, outputPS] = ...                 
                patternsearch(@(p_est) sum(reshape(myModel(p_est*p.D), ...
                p.I*p.J, 1).^2), ...
                p_est, [], [], [], [], p_lbPS, p_ubPS , [], optionsPS);
            calcTimePS = toc(startTime) - calcTimeLSQ;
            
            [~, s.SSE, s.Residuals, ~, ~, ~, Jacobian] = ...
                lsqnonlin(@(p_est) myModel(p_est*p.D), ...
                p_est, p_lb, p_ub, optimoptions(optionsLSQ, 'Display', ...
                'off', 'MaxFunctionEvaluations', 0));            
        case 0
            exitflagPS = [];
            outputPS = [];
            calcTimePS = [];
    end
    
    p_est = p_est*p.D;
    p.HGcalc = n18rvcalc./(n18rvcalc + n0rvcalc);
    
    %% statistics
    s.calcTimeLSQ = calcTimeLSQ;                                            % calculation time for lsqnonlin
    s.calcTimePS = calcTimePS;                                              % calculation time for patternsearch
    s.exitflagLSQ = exitflagLSQ;
    s.exitflagPS = exitflagPS;
    s.outputLSQ = outputLSQ;
    s.outputPS = outputPS;
    s.Nobs = p.I*p.J;                                                       % number of observations/data points (coefficients)
    s.NumCoefficients = length(p_est);                                      % number of coefficients
    s.CoefficientNames = CoefficientNames;
    s.VariableNames = {'t' 'n18' 'n0' 'Tm' 'Tr' 'mc' 'pr' ...
        'dn18v' 'dn0v' 'n18rv'}';
    s.NumVariables = length(s.VariableNames);                               % number of variables  
    
    % help: edit nlinfit (ll. 364 - 378 & 390 - 421)
    %       edit nonlinearmodel (ll. 1145 & 1165)    
    s.J = full(Jacobian);                                                   % Jacobi-matrix (coefficients)    
    [~, R] = qr(s.J, 0);                                                    % QR factorization of Jacobi matrix -> s.J = Q*R
    if s.Nobs <= s.NumCoefficients || condest(R) > 1/sqrt(eps(class(p_est)))
        if s.Nobs <= s.NumCoefficients
            warning(message('stats:nlinfit:Overparameterized'));
        elseif condest(R) > 1/sqrt(eps(class(p_est)))
            if any(all(abs(s.J)<sqrt(eps(norm(s.J,1))),1),2)                % one or more columns of zeros
                warning(message('stats:nlinfit:ModelConstantWRTParam'));
            else                                                            % no columns of zeros -> estimates highly correlated        
                warning(message('stats:nlinfit:IllConditionedJacobian'));
            end
        end
        TolSVD = eps(class(p_est));
        [~, rankJ, pinvJTJ] = internal.stats.isEstimable(eye(numel(p_est)),...
            'DesignMatrix', s.J, 'TolSVD', TolSVD);
    else
        rankJ = size(s.J,2);
        Rinv = R\eye(size(R));                                              % inverse of R
        pinvJTJ = Rinv*Rinv';                                               % pseudo inverse of J'*J.
    end
    s.DFE = s.Nobs - rankJ;                                                 % error degree of freedom (coefficients)
    s.NumEstimatedCoefficients = s.Nobs - s.DFE;                            % number of estimated coefficients
    
    s.MSE = s.SSE/s.DFE;                                                    % mean squared error (coefficients)
    s.RMSE = sqrt(s.MSE);                                                   % root mean squared error (coefficients)
    s.CoefficientCovariance = s.MSE*pinvJTJ;                                % variance - covariance matrix for coefficient estimates (coefficients)

    % help: edit nonlinearmodel (ll. 1183 - 1211)
    n18rvmean = sum(p.weight.*mb.n18rv)/sum(p.weight);                      % weighted mean of mb.n18rv
    s.SST = sum(p.weight(:).*(mb.n18rv(:) - n18rvmean).^2);                 % total sum of squares (coefficients)
    s.SST0 = sum(p.weight(:).*(mb.n18rv(:) - 0).^2);
    s.SSR = sum(p.weight(:).*(n18rvcalc(:) - n18rvmean).^2);                % regression sum of squares (coefficients)

%     n18rvmean = mean(mb.n18rv, 'all');                                      % mean of mb.n18rv
%     s.SST = sum((mb.n18rv - n18rvmean).^2, 'all');                          % total sum of squares (coefficients)
%     s.SST0 = sum(p.weight.*(mb.n18rv - 0).^2, 'all');
%     s.SSR = sum(p.weight.*(n18rvcalc - n18rvmean).^2, 'all');               % regression sum of squares (coefficients)
       
    % help: edit classreg.regr.modelutils.rsquared (ll. 28 - 39)
    %       edit nonlinearmodel (ll. 1234 - 1251)
    s.DFT = s.Nobs - 1;                                                     % corrected degrees of freedom (coefficients = total)
    s.Rsquared.Ordinary = 1-(s.SSE./s.SST);                                 % Ordinary (unadjusted) R-squared (model)
    s.Rsquared.Adjusted =  1-(s.SSE./s.SST)*(s.DFT./s.DFE);                 % Adjusted R-squared, adjusted for the number of coefficients (model)
    s.LogLikelihood = -(s.DFE + s.Nobs*log(2*pi) + ...                      % log of likelihood function at coefficient estimates (coefficients)
        sum(log((s.DFE/s.Nobs*s.MSE)./p.weight(:))))/2;                     % the higher the value the better the regression 
   
    % help: edit nonlinearmodel (ll. 719 - 724)
    s.alpha = 0.05;                                                         % alpha (coefficients)
    s.SE = p.D*sqrt(diag(s.CoefficientCovariance));                         % standard error of the estimate (coefficients)
    s.delta = s.SE*tinv(1-s.alpha/2, s.DFE);                                % half width of confidence interval (coefficients)
    s.ci = [p_est' - s.delta p_est' + s.delta];                             % (1-alpha) confidence interval (coefficients)

    % help: edit nonlinearmodel (ll. 805 - 834 & 1214 - 1224)
    s.Junw = bsxfun(@ldivide, sqrt(p.weight(:)), s.J);                      % unweighted Jacobi-matrix (coefficients)
    Jmin = min(s.Junw,[],1);
    Jmax = max(s.Junw,[],1);
    if any( abs(Jmax-Jmin) <= ...                                           % constant column of s.Junw
            sqrt(eps(class(s.Junw))) * (abs(Jmax) + abs(Jmin)) ) ...
            &&  s.NumEstimatedCoefficients > 1
        % Compare full model vs. intercept only model
        s.ssr = s.SST - s.SSE;
        s.dfr = s.NumEstimatedCoefficients - 1;
        s.dfe = s.Nobs - 1 - s.dfr;
    else
        % Compare full model vs. zero model
        s.ssr = s.SST0 - s.SSE;                                             % regression sum of squares (model)
        s.dfr = s.NumEstimatedCoefficients;                                 % corrected degrees of freedom (model)                   
        s.dfe = s.Nobs - s.dfr;                                             % error degree of freedom (model)
    end
    s.f = (s.ssr./s.dfr) / (s.SSE/s.dfe);                                   % F statistic for a test that the model is intercept only / zero (model)
    s.p = fcdf(1./s.f, s.dfe, s.dfr);                                       % p-value for the F statistic (model)

    % help: edit nonlinearmodel (ll. 841 - 842)
    %       edit classreg.regr.modelutils.tstats (ll. 43 - 45)
    s.tStat = p_est'./s.SE;                                                 % t statistic for a test that the coefficient is zero (coefficients)                      
    s.pValue = 2*tcdf(-abs(s.tStat), s.DFE);                                % p-value for the t statistic (coefficients)
    
    % significance and results of null hypothesis test (coefficients)
    s.signif = cell(s.NumCoefficients,1);
    s.NonZero = cell(s.NumCoefficients,1);
    for k = 1:s.NumCoefficients
        % Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
        if 0 <= s.pValue(k) && s.pValue(k) <= 0.001
            s.signif(k) = cellstr(char('***'));
        elseif 0.001 < s.pValue(k) && s.pValue(k) <= 0.01
            s.signif(k) = cellstr(char('**'));
        elseif 0.01 < s.pValue(k) && s.pValue(k) <= 0.05
            s.signif(k) = cellstr(char('*'));
        elseif 0.05 < s.pValue(k) && s.pValue(k) <= 0.1
            s.signif(k) = cellstr(char('.'));
        elseif s.pValue(k) > 0.1
            s.signif(k) = cellstr(char(''));
        end
        if s.pValue(k) <= s.alpha
            s.NonZero(k) = cellstr(char('yes'));
        else
            s.NonZero(k) = cellstr(char(''));
        end
    end     
        
    if isempty(char(s.signif)) || isempty(char(s.NonZero))
        s.Coefficients = table(p_est', s.SE, s.tStat, s.pValue, ...
            'RowNames', CoefficientNames, 'VariableNames',...
            {'Estimate', 'SE', 'tStat','pValue'});
    else
        s.Coefficients = table(p_est', s.SE, s.tStat, s.pValue, ...
            char(s.signif), char(s.NonZero), ...
            'RowNames', CoefficientNames, 'VariableNames',...
            {'Estimate', 'SE', 'tStat','pValue', 'signif', 'NonZero'});
    end
    
    % diagnostics
    % https://de.mathworks.com/help/stats/nonlinearmodel.plotdiagnostics.html#btd2jcc-1
    % help: edit nonlinearmodel (ll. 909 - 936)
    [Q, ~, ~] = qr(s.J, 0);
    d.HatMatrix = ...
        bsxfun(@times, bsxfun(@times, 1./sqrt(p.weight(:)), ...
        Q(:,1:s.NumEstimatedCoefficients)*Q(:,1:s.NumEstimatedCoefficients)'), ...
        sqrt(p.weight(:))');

    % help: edit nonlinearmodel (ll. 1157 - 1180 & 938 - 947)
    d.Leverage = sum(abs(Q(:,1:s.NumEstimatedCoefficients)).^2, 2);
    d.CooksDistance = s.Residuals(:).^2./...                                % Cook's distance is useful for identifying outliers in the experimental values
        (s.NumCoefficients*s.MSE).*(d.Leverage(:)./(1-d.Leverage(:)).^2);
    d.CooksDistanceIndex = find((d.CooksDistance)>3*mean(d.CooksDistance)); % observations with Cooks distance values that exceed the recommended threshold value   
        
    %% save the estimated parameters in workspace    
    % model
    mdl.MSE = s.MSE;
    mdl.Residuals = s.Residuals./sqrt(p.weight);    
    mdl.Fitted = n18rvcalc;
    mdl.RMSE = s.RMSE;
    mdl.Diagnostics = table(d.Leverage, d.CooksDistance, d.HatMatrix, ...
        'VariableNames', {'Leverage', 'CooksDistance', 'HatMarix'});
    mdl.WeightedResiduals = s.Residuals;
    mdl.NumVariables = s.NumVariables;
    mdl.VariableNames = s.VariableNames;
    mdl.NumPredictors = mdl.NumVariables - 1;
    mdl.PredictorNames = mdl.VariableNames(1:mdl.NumPredictors,1);
    mdl.ResponseName = mdl.VariableNames(end,1);
    mdl.NumObservations = s.Nobs;
    mdl.ObservationInfo = table(p.weight(:), zeros(s.Nobs,1), ...
        zeros(s.Nobs,1), ones(s.Nobs,1), ...
        'VariableNames', {'Weights', 'Excluded', 'Missing', 'Subset'});
    mdl.Variables = table(repmat(p.t,(s.Nobs./p.I),1), ...
        mb.n18(:), mb.n0(:) , ...        
        repmat(p.Tm(1,1),s.Nobs,1), repelem(p.Tr(:),p.I,1), ...
        reshape(repmat(p.mc,p.I,1),[s.Nobs,1]), repelem(p.pr(:),p.I,1),  ...
        mb.dn18v(:), mb.dn0v(:), ...
        mb.n18rv(:), ...
        'VariableNames', s.VariableNames');
    mdl.ObservationNames = {};
    mdl.Formula.Expression = Rate;    
    mdl.LogLikelihood = s.LogLikelihood;
    mdl.DFE = s.DFE;
    mdl.SSE = s.SSE;
    mdl.SST = s.SST;
    mdl.SSR = s.SSR;
    mdl.CoefficientCovariance = s.CoefficientCovariance;
    mdl.CoefficientNames = s.CoefficientNames;
    mdl.NumCoefficients  = s.NumCoefficients;
    mdl.NumEstimatedCoefficients  = s.NumEstimatedCoefficients;
    mdl.Coefficients = s.Coefficients;
    mdl.Rsquared.Ordinary = s.Rsquared.Ordinary;
    mdl.Rsquared.Adjusted = s.Rsquared.Adjusted;
    
    %% reaction paramters
    r.EA  = p_est(cell2sym(CoefficientNames) == 'A').*...                   % activation energy (J mol^-1)
        (const.R.*p.Tm(1,1));      
    r.k0 = p_est(cell2sym(CoefficientNames) == 'k').*...                    % frequency of collision (mol^(1-n) s^-1  l^(n-1) frate^-n), for unit of frate see myReactionRate.m
        (exp(p_est(cell2sym(CoefficientNames) == 'A')));
    
    %% plot results
    ymin = 0;
    ymax = 1;
    legendname = {'experiment', 'model'};
    legendnameP = {'1 bar', '2 bar', '3 bar', '4 bar', '5 bar'};
    legendnameT = {'270 °C','280 °C', '290 °C', '300 °C', '310 °C'};
    xlabelname = 't';
    xlabelunit = '(s)';
    ylabelname = 'HG';
    ylabelunit = '(-)';

    % HG vs HGcalc
    figure()
    hold on
    Colors = get(gca, 'ColorOrder');
    Markers = {'s','d','o','.','x'};
    h = zeros(p.nP, p.nTheta);
    for j = 1:p.nP
        dj = linspace(p.nTheta*(j-1)+1, j*p.nTheta, p.nTheta);
        if j == 1
            line([ymin ymax],[ymin ymax], 'Color', 'k');
            
            xlim([ymin,ymax])
            ylim([ymin,ymax])
            
            xlabel(strcat(ylabelname, ' (exp)', {' '} , ylabelunit))   
            ylabel(strcat(ylabelname, ' (calc)', {' '} , ylabelunit))  
        end
        h(j,:) = plot(p.HG(:,dj), p.HGcalc(:,dj),'x', 'Color', Colors(j,:));        
        for zz = 1:p.nTheta
            set(h(j,zz),'Marker', Markers{zz})
        end           
    end
    if p.nP == p.nTheta
        h = diag(h);        
        legendnameP = legendnameP(p.pr(1:p.nTheta:p.J));
        legendnameT = legendnameT(5:-1:(6-p.nTheta));
    elseif p.nP < p.nTheta
        if p.nP ~= 1
            h = [diag(h); h(p.nP, p.nP+1:p.nTheta)'];
        end
        legendnameP = [legendnameP(p.pr(1:p.nTheta:p.J)), ...
            cell(1, p.nTheta-p.nP)];
        legendnameP(p.nP+1:p.nTheta) = legendnameP(p.nP);
        legendnameT = legendnameT(5:-1:(6-p.nTheta));
    elseif p.nP > p.nTheta
        if p.nTheta ~= 1
            h = [diag(h); h(p.nTheta+1:p.nP,p.nTheta)'];
        end
        legendnameP = legendnameP(p.pr(1:p.nTheta:p.J));
        legendnameT = [legendnameT(5:-1:(6-p.nTheta)), cell(1, p.nP-p.nTheta)];
        legendnameT(p.nTheta+1:p.nP) = legendnameT(p.nTheta);
    end
    [l, ~, plots] = legend(gca, h, strcat(legendnameP, {' & '}, ...
        legendnameT), 'Location', 'northwest');
    for idx = 1:length(l.String)
        if idx <= p.nP
            l.String{idx} = ...
                ['{\color[rgb]{' num2str(plots(idx).Color) '} ' ...
                l.String{idx}(1:5) '}' l.String{idx}(6:end)];
        else
            l.String{idx} = ...
                ['{\color[rgb]{' num2str([1 1 1]) '} ' ...
                l.String{idx}(1:7) '}' l.String{idx}(8:end)];
        end
    end
    hold off
    
    if p.plot == 1
        figure('units','normalized','outerposition',[0 0 1 1])
        for j = 1:p.nP
            dj = linspace(p.nTheta*(j-1)+1, j*p.nTheta, p.nTheta);
            if p.nP > 1
                subplot(2,round(p.nP/2),j)
            end
            hold on
            h1 = plot(p.t, p.HGcalc(:,dj),'-');
%             ax = gca;
%             ax.ColorOrderIndex = 1;
            set(gca,'ColorOrderIndex',1);
            h2 = plot(p.t, p.HG(:,dj),'x');
            hold off
            title(sprintf('\\bf %d bar',p.pr((j-1)*p.nTheta+1)))
            ylim([ymin,ymax])

            if j == 1
                legend(h1(p.nTheta:-1:1), legendnameT(end:-1:1))
                ylabel(strcat(ylabelname, {' '} , ylabelunit)) 
                if p.nP == 1
                    xlabel(strcat(xlabelname, {' '} , xlabelunit))                
                end
            elseif j == 2
                legend([h2(p.nTheta) h1(p.nTheta)], legendname)
            elseif j == 2 && p.nP <= 3
                xlabel(strcat(xlabelname, {' '} , xlabelunit))          
                if p.nP <= 2
                    ylabel(strcat(ylabelname, {' '} , ylabelunit))                  
                end
            elseif j == 3
                xlabel(strcat(xlabelname, {' '} , xlabelunit))
                if p.nP <= 4
                     ylabel(strcat(ylabelname, {' '} , ylabelunit))                
                end
            elseif j == 4
                xlabel(strcat(xlabelname, {' '} , xlabelunit))
                if p.nP >= 5
                    ylabel(strcat(ylabelname, {' '} , ylabelunit)) 
                end
            elseif j == 5
                xlabel(strcat(xlabelname, {' '} , xlabelunit))
            end
        end

        % plotDiagnostics(mdl, 'cookd')
        figure()
%         ax = gca;
%         ax.XGrid = 'on';
       set(gca,'XGrid','on');
%         ax.XTick = 0:p.I:mdl.NumObservations;
       set(gca,'XTick',0:p.I:mdl.NumObservations);
        hold on
        plot(1:mdl.NumObservations, d.CooksDistance(:), 'rx')
        plot(1:mdl.NumObservations, ...
            3*mean(d.CooksDistance(:))*ones(mdl.NumObservations,1), 'k--')
        strValues = strtrim(cellstr(num2str(d.CooksDistanceIndex(:),'%d')));
        text(d.CooksDistanceIndex(:), d.CooksDistance(d.CooksDistanceIndex(:)), ...
            strValues,'VerticalAlignment','top', 'HorizontalAlignment','left');
        for j = 1:p.nP
%            line([p.I*p.nTheta*j p.I*p.nTheta*j], ...
%                [0 round(max(d.CooksDistance(:)), ...
%                -floor(log10(max(d.CooksDistance(:)))))], 'Color', 'red');
            line([p.I*p.nTheta*j p.I*p.nTheta*j], ylim, 'Color', 'red');
        end
        hold off        
        title('Case order plot of Cook''s distance')
        xlabel('Row number')
        ylabel('Cook''s distance') 

        % plotResiduals(mdl, 'probability')
        figure()
        probplot(mdl.Residuals(:))
    end
    
    %% display diagnostics & calculation time
    fprintf('Nonlinear regression model: \n \t rate ~ %s \n', ...
        char(mdl.Formula.Expression))
    disp(s.Coefficients);
    fprintf('Number of observations: %d, Error degrees of freedom: %d \n', ...
        s.Nobs, s.DFE)
    fprintf('Root Mean Squared Error: %d \n', s.RMSE)
    fprintf('R-Squared: %d, Adjusted R-Squared: %d \n', ...
        s.Rsquared.Ordinary, s.Rsquared.Adjusted)
    if s.dfr == s.NumEstimatedCoefficients
        fprintf('F-statistics vs. zero model: %d, p-value = %d \n', s.f, s.p)
    else
        fprintf('F-statistics vs. intercept only model: %d, p-value = %d \n', ...
            s.f, s.p)
    end
    
    %% nested functions
    function Residuals = myModel(p_est)
        J = 2*p.J;
        JPattern = spdiags([ones(J,1) ones(J,1) ones(J,1)], [-J/2, 0, J/2], J, J);
        options = odeset('NonNegative', 1:1:J, 'JPattern', JPattern);

        [~, nrvcalc] = ode15s(@(t, nrv, p_est) myMolarChange(t, nrv, ...
            p_est), p.t', [mb.n18(1,:) mb.n0(1,:)], options, p_est);

        if isfield(p,'warnMH')
            p = rmfield(p,'warnMH');
        end
        
        n18rvcalc = nrvcalc(:,1:p.J);
        n0rvcalc = nrvcalc(:,p.J+1:end); 
        Residuals = (mb.n18rv - n18rvcalc).*sqrt(p.weight);
        
        function dndt = myMolarChange(~, nrv, p_est)            
            n18 = nrv(1:p.J)';                                              % molar amount of H18-DBT in liquid phase (mol)
            n0 = nrv(p.J+1:end)';                                           % molar amount of H18-DBT in liquid phase (mol)  
            p.p_est = p_est;
            
            [reff4, V, ~, ~, p] = myReactionRate(const, p, p.pr, ...        % effective reaction rate (mol l_bed^-1 s^-1) & volume of liquid phase (l)
                p.Tr, n18, n0);     
            dn18dt = (const.nu(1).* ...                                     % molar change of H18-DBT in liquid phase (mol s^-1)
                (V+p.Vc.*(1-const.epsBed)/1000).*reff4)';
            dn0dt = (const.nu(2).*...                                       % molar change of H0-DBT in liquid phase (mol s^-1)
                (V+p.Vc.*(1-const.epsBed)/1000).*reff4)'; 

%             [reff4, ~, ~, ~, p] = myReactionRate(const, p, p.pr, ...        % effective reaction rate (mol l_bed^-1 s^-1) & volume of liquid phase (l)
%                 p.Tr, n18, n0);     
%             dn18dt = (const.nu(1).*(p.Vc/1000).*reff4)';                    % molar change of H18-DBT in liquid phase (mol s^-1)
%             dn0dt = (const.nu(2).*(p.Vc/1000).*reff4)';                     % molar change of H0-DBT in liquid phase (mol s^-1)

            dndt = [dn18dt; dn0dt];
        end
    end
end