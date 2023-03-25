function [x, Tr, HG, eta, etaexp, residual, exitflag]=etaPPDSFit
    HG = [0 4.88 11.81 19.55 26.54 29.87 37.58 43.78 49.08 54.82 63.65 ...
        68.38 77.30 88.43 96.35 100]'./100;
    Tr = [20 30 40 50 60 70	80 100 120 140] + 273.15;
    etaexp=[51.09	28.25	17.32	11.64	8.30	6.21	4.81	3.16	2.24	1.68
            53.51	29.33	17.90	11.95	8.48	6.31	4.88	3.19	2.26	1.69
            60.70	32.69	19.60	12.90	9.05	6.69	5.13	3.32	2.33	1.74
            72.42	37.73	22.15	14.34	9.92	7.24	5.50	3.49	2.43	1.80
            95.73	47.61	26.73	16.76	11.32	8.11	6.07	3.78	2.59	1.89
            93.94	46.78	26.36	16.61	11.21	8.05	6.02	3.74	2.56	1.88
            135.00	63.61	34.07	20.60	13.49	9.41	6.90	4.17	2.79	2.01
            191.10	85.09	43.61	25.48	16.08	10.94	7.84	4.59	3.01	2.13
            268.30	114.20	55.81	31.25	19.17	12.72	8.94	5.08	3.26	2.28
            268.70	113.40	55.33	30.98	19.04	12.60	8.85	5.02	3.22	2.25
            631.50	242.40	106.10	53.97	30.48	18.86	12.50	6.51	3.96	2.65
            407.30	162.40	74.95	40.05	23.64	15.19	10.40	5.67	3.54	2.42
            408.10	160.70	73.62	39.22	23.17	14.90	10.20	5.60	3.50	2.40
            470.10	180.40	80.87	42.43	24.77	15.76	10.70	5.82	3.61	2.47
            494.10	184.30	82.01	42.93	25.02	15.93	10.81	5.84	3.63	2.47
            413.20	157.90	71.56	38.10	22.55	14.54	10.00	5.52	3.47	2.39];
    
    Tr = Tr(3:end);
    HG = HG(1:end);
    etaexp=etaexp(1:end,3:end);
    
    eta18expTVT = exp(-430.448 + 25938.81./Tr + 61.220.*log(Tr));
    eta12expTVT = exp(-463.825 + 29039.182./Tr + 65.505.*log(Tr));
    eta6expTVT = exp(-441.778 + 26446.189./Tr + 62.888.*log(Tr));
    eta0expTVT = exp(-275.593 + 16793.053./Tr + 39.122.*log(Tr)); 
    etaexpTVT = [eta18expTVT; eta12expTVT; eta6expTVT; eta0expTVT];
    
    A = []; B = []; C = []; D = []; E = [];
    N = 6;
    optionsLSQ = optimoptions('lsqnonlin', ...
        'Display', 'iter', 'FiniteDifferenceType', 'central', ...
        'MaxFunctionEvaluations', 100000, 'MaxIter', 100000, ...
        'TolFun', 1e-10);
    x0 = [-31, -25, 250, 250, 1e-8];
    lb = [-100, -100, -1000, -1000, 0];
    ub = [+100, +100, 1000, +1000, 1e-2];
    x0 = reshape([x0', zeros(size(x0,2),N-1)]',1,size(x0,2)*N);
    lb = reshape([lb', -Inf(size(lb,2),N-1)]',1,size(lb,2)*N);
    ub = reshape([ub', Inf(size(ub,2),N-1)]',1,size(ub,2)*N);
    [x, resnorm, residual, exitflag] = lsqnonlin(@(x) fun(x), x0, lb, ub, optionsLSQ);
%     x = reshape(x, N, 5);

    eta = E.*exp(A.*-((Tr - C)./(Tr - D)).^(1/3) + ...
            B.*-((Tr - C)./(Tr - D)).^(1/3).*((C - Tr)./(Tr - D)));

%     ydash = mean(etaexp(:));
%     R2 = sum((eta(:)-ydash).^2)./sum((etaexp(:)-ydash).^2);
    R2 = 1 - resnorm./sum((etaexp(:)-mean(etaexp(:))).^2);
    R2Adj = 1-(1-R2).*(numel(etaexp)-1)./(numel(etaexp)-30);

    hold on
    Colors = (1-((1:size(etaexp,1))./size(etaexp,1)))'.*ones(1,3);
    h = zeros(size(etaexp,1),1);
    hexp = zeros(size(h));
    for i = 1:size(etaexp,1)    
        h(i) = plot(Tr, eta(i,:),'Color', Colors(i,:), 'LineStyle', '-');    
        hexp(i) = plot(Tr, etaexp(i,:), 'Color', Colors(i,:), ...
            'Marker', 'x', 'LineStyle', 'none');
        if i == size(etaexp,1)
            hexpTVT = plot(Tr, etaexpTVT, 'gx-');
        end
    end
    legend([h;hexpTVT(1)] , [num2str(round(100*HG,1)); ' TVT'])
    
    for i = 1:numel(HG)
        figure() 
        plot(Tr, eta(i,:),'b-', Tr, etaexp(i,:),'rx')
    end
    
    function y = fun(x)
        A = modelfun(x(1:N));
        B = modelfun(x(N+1:2*N));
        C = modelfun(x(2*N+1:3*N));
        D = modelfun(x(3*N+1:4*N));
        E = modelfun(x(4*N+1:5*N));
        etaPPDS = E.*exp(A.*-((Tr - C)./(Tr - D)).^(1/3) + ...
            B.*-((Tr - C)./(Tr - D)).^(1/3).*((C - Tr)./(Tr - D)));
        y = etaPPDS - etaexp;
        function y = modelfun(x)
            y = x(1) + x(2).*HG + x(3).*HG.^2 + x(4).*HG.^3 + ...
                x(5).*mod(HG, 2/3) + x(6).*mod(HG, 2/3).^2;
        end
    end
end