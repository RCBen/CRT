function SSE = error_2xz_400g
    mydir = pwd;
    idcs = strfind(mydir,'\');
    p.subdir = mydir(1:idcs(end-3)-1);

    % constants
    if isempty(which('getConstants'))      
        addpath(strcat(p.subdir,'\functions'))
    end
    const = getConstants;


    %% experimental data
    name = {'D102' 'D101/103' 'D105' 'D114'};
    exp = {'W' 'Z' 'W' 'W'};
    mc = 400*ones(1,17);
    Q = [9.99 15.02 10 10];
    pr = [2 1.49 1.05 2.01];
    Tr = [279.65 293.98 278.91 307.64]+273.15;  
    
    % F6 to F10 from 200g DoE experiments
    F6 = [283.96 300.81 284.55 315.33];
    F7 = [286.41 303.30 286.90 318.04];
    F8 = [286.34 304.07 287.43 318.19];
    F9 = [286.20 304.72 287.32 317.87];
    F10 = [285.44 306.16 286.08 318.24];    

%   data based QH2 measurements via MFM
%     HGout = ones(size(name));
    QH2out = [1373 3378 1703 3214];                                         % experimental values of module 2
    
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
    
    QPlan = [10 15 10 10];
    prPlan = [2 1.5 1 2];
    TrPlan = [280 295 280 310];
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

%     title('Experiment vs. simulation');
% %     xlabel(name);
%     ylabel('hydrogenation grade [%]');
%     ylim([0 100])
%     
%     bar([HGoutcalc'*100 HGout'*100])
%     legend('calc', 'exp');
    
    title('Experiment vs. simulation');
    ylabel('hydrogen flow rate [Nl/min]');
    ylim([0 5])
    
    bar([QH2outcalc'/1000 QH2out'/1000])
    legend('calc', 'exp');
end