function SSE = error_2xz
    mydir = pwd;
    idcs = strfind(mydir,'\');
    p.subdir = mydir(1:idcs(end-3)-1);

    % constants
    if isempty(which('getConstants'))      
        addpath(strcat(p.subdir,'\functions'))
    end
    const = getConstants;


    %% experimental data
    name = {'D101' 'D102' 'D103' 'D104' 'D105' 'D106' 'D108' 'D109' ...
        'D113' 'D114' 'D115' 'D116' 'D117' 'X101.1' 'X101.2' 'X102.1' ...
        'X102.2'};
    exp = {'Z' 'W' 'Z' 'W' 'W' 'W' 'W' 'W' 'Z' 'W' '*' '*' 'W' ...
        'Val' 'Val' 'Val' 'Val'};
    mc = 200*ones(1,17);
    Q = [15.02 9.97 14.97 19.98 9.99 10.00 19.98 19.99 14.99 9.98 ...
        14.99 15.02 20.01 10.31 15.53 18.01 11.54];
    pr = [1.52 1.99 1.49 1.04 1.05 1.06 1.95 1.97 1.47 1.95 2.23 1.04 ...
        1.05 1.94 1.79 1.12 1.53];
    Tr = [295.31 279.69 295.45 279.76 279.93 310.22 308.47 279.03 ...
        294.34 308.82 295.04 294.03 307.30 284.45 287.06 293.65 301.44] + ...
        273.15;    
    
    F6 = [300.91 283.96 300.81 285.85 284.55 315.68 317.85 285.54 ...
        300.81 315.33 301.13 301.95 317.05 289.53 294.46 302.31 309.25];
    F7 = [303.29 286.41 303.30 287.37 286.90 319.14 319.88 286.43 ...
        302.77 318.04 303.77 302.93 318.86 291.71 295.47 302.62 310.82];
    F8 = [304.00 286.34 304.07 287.55 287.43 319.48 320.09 285.38 ...
        302.86 318.19 304.29 302.18 318.28 291.40 294.09 301.98 310.62];
    F9 = [304.72 286.20 304.72 288.62 287.32 318.28 320.66 286.17 ...
        303.64 317.87 304.91 302.66 317.96 291.36 294.12 303.04 310.62];
    F10 = [306.01 285.44 306.16 288.97 286.08 317.96 321.11 288.43 ...
        305.53 318.24 306.11 304.80 318.77 291.63 296.95 306.42 311.83];    

%   data based QH2 measurements via MFM
% %     HGout = [0.7418 0.8595 0.7383 0.8417 0.7989 0.6249 0.6820 0.9015 ...
% %         0.7680 0.6100 0.7906 0.7520 0.6924 0.8522 0.8544 0.7779 0.6643];
    QH2out = [2642 954 2673 2142 1356 2538 4335 1315 2370 2653 2138 ...     % experimental values of module 2
        2538 4198 1050 1560 2703 2626];
    
%   data based on liquid samples and NMR measurements
% %     HGout = [0.7757 0.8734 0.7720 0.8661 0.8282 0.6788 0.7148 0.8605 ...
% %         0.8015 0.6502 0.8199 0.7895 0.7286 0.8402 0.8716 0.8116 0.7194];
%     QH2out = [2295 860 2329 1810 1156 2170 3888 1874 2028 2380 1838 ...     % experimental values of module 2
%         2154 3703 1135 1377 2289 2192];
%     QH2out2 = 1.15*(2505.1879 + 608.1759.*(Q-15)/5 - ...                    % calculated, experimental values with DoE equation of module 1 (times 1.15 (15 percent improvement))
%         140.0390.*(pr-1.5)./0.5 + 1027.2766.*(Tr-(295+273.15))/15 + ...
%         302.9439.*(Q-15)/5.*(Tr-(295+273.15))/15 + ...
%         188.1855.*(pr-1.5)/0.5.*(Tr-(295+273.15))/15);

    HGout = 1-QH2out./Q./(getRho(const, pr, 293.15, 1)/1000)*const.M18/...
        const.nu(3)/const.MH2*const.rhoH2n/1000;
       
    HGoutcalc = zeros(size(name));
    QH2outcalc = zeros(size(name));
    Trcalc = zeros(size(name));
    nx = zeros(size(name));
    nz = zeros(size(name));
% 
%     f = waitbar(0,'1','Name','progress bar...',...
%         'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    for i = 1:numel(name)
        Tw = [F6(i) mean([F6(i), F7(i)]) F7(i) mean([F7(i), F8(i)]) F8(i) ...
            mean([F8(i), F9(i)]) F9(i) mean([F9(i), F10(i)]) F10(i)] + ...
            273.15;
        [HGoutcalc(i), QH2outcalc(i), prcalc, Trcalc(i)] = ...
            main_2Dxz(mc(i), Q(i), pr(i), Tr(i), Tw);

        nx(i) = size(prcalc,2);
        nz(i) = size(prcalc,1);

%         % Check for clicked Cancel button
%         if getappdata(f,'canceling')
%             break
%         end
%         
%         % Update waitbar and message
%         waitbar(i/numel(name),f,sprintf('%d of %d finished', i, numel(name)))
    end
%     delete(f)
    
    
    QPlan = [15 10 15 20 10 10 20 20 15 10 15 15 20 10.3 15.5 18 11.5];
    prPlan = [1.5 2 1.5 1 1 1 2 2 1.5 2 2.21 1 1 1.96 1.76 1.14 1.53];
    TrPlan = [295 280 295 280 280 310 310 280 295 310 295 295 310 ...
        284.7 288.4 294.6 302.7];
    ErrorHG = (HGoutcalc - HGout)./HGout;
    ErrorQH2 = (QH2outcalc - QH2out)./QH2out;
    dTr = (Trcalc - Tr)';
    SSE = sum(ErrorQH2(:).^2);
    
    tbl = table(QPlan', prPlan', TrPlan', dTr, HGout', HGoutcalc', ...
        ErrorHG'*100, QH2out', ceil(QH2outcalc'),ErrorQH2'*100, ...
        'RowNames', strcat(name, {' ('}, exp, {')'}), ...
        'VariableNames', {'VLOHC', 'pr', 'Tr', 'Trcalc-Tr', 'HGout', ...
        'HGoutcalc', 'ErrorHG', 'QH2out', 'QH2outcalc', 'ErrorQH2'});
    disp(tbl)

    figure()
    grid on
    hold on 

    title('Experiment vs. simulation');
%     xlabel(name);
    ylabel('hydrogenation grade [%]');
    ylim([0 100])
    
    bar([HGoutcalc'*100 HGout'*100])
    legend('calc', 'exp');
    
%     title('Experiment vs. simulation');
%     ylabel('hydrogen flow rate [Nl/min]');
%     ylim([0 5])
%     
%     bar([QH2outcalc'/1000 QH2out'/1000])
%     legend('calc', 'exp');
end