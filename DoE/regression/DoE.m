function [mdl, Error] = DoE(module, result)
    % [mdl, Error] = DoE(module, result)
    %
    % module            1   PDBoden V1
    %                   2   PDBoden V2
    % result            1   HG      (hydrogenation grade)
    %                   2   QH2     (hydrogen volume flow rate)
    %                   3   PthH2   (thermal power output due to QH2)
    %                   4   Etar    (energetic efficiency of dehydrogenation)
    %                   5   Eta     (energetic efficiency wo pre-heating)
    %                   6   Etaw0   (energetic efficiency with pre-heating)
    %                   7   LB      (change in light boiling products)
    %                   8   HB      (change in high boiling products)
    %                   9   RTD     (residence time distribution)

    if nargin == 0
        module = 1;
        result = 1;
    end
    plotting = 0;

    %% DoE factor levels
    % factor levels
    Qmin = 10;
    Qmax = 20;
    pmin = 1; 
    pmax = 2; 
    Tmin = 280;
    Tmax = 310;
    
    % normalized factor levels
    nmin = -1;
    nmax = +1;
    
    % Linear equation (slope, intercept)
    mQ = (nmax-nmin)/(Qmax-Qmin);
    tQ = -mQ*mean([Qmax, Qmin]);
    
    mp = (nmax-nmin)/(pmax-pmin);
    tp = -mp*mean([pmax, pmin]);
    
    mT = (nmax-nmin)/(Tmax-Tmin);
    tT = -mT*mean([Tmax, Tmin]);
    
    %% Read in data
    file = '2018 DoE Ergebnisse.xlsx';
    sheet = 'Ergebnis';
    Data = xlsread(fullfile(pwd, 'Messdaten', file), sheet, 'C3:AI38');
    DataVal = xlsread(fullfile(pwd, 'Messdaten', file), sheet, 'C45:AI52');
    
    %% Process data
    Variables = {'QLOHC' 'p' 'T' 'HG' 'QH2' 'PthH2' 'Etar' 'Eta' 'Etaw0' ...
        'LB' 'HB' 'RTD'};
    Units = {'[ml min^-1]' '[bara]' '[°C]' '[-]' '[Nml min^-1]' ...
        '[W]' '[-]' '[-]' '[-]' '[-]' '[-]' '[s]'};
    yRange = [1:3, 8, 12, 19, 25:27, 31:33];
    
    Factors = Variables(1:3);
    Response = Variables(result+3);
    Unit = Units(result+3);
    
    switch module
        case 1
            modulename = 'M1';
            xRange = 1:18;
            xRangeVal = 1:4;
        case 2
            modulename = 'M2';
            xRange = 18+[1:6, 8:9, 13:17];
            xRangeVal = 5:8;
    end
    Data = Data(xRange,yRange);
    Qn = polyval([mQ tQ], Data(:,1));
    pn = polyval([mp tp], Data(:,2));
    Tn = polyval([mT tT], Data(:,3));
    Result = Data(:, result+3);

    %% Estimate the coefficients of the model
    VariableNames = [Factors Response];
    tbl = table(Qn, pn, Tn, Result, 'VariableNames', ...
        VariableNames);
    
    mdl = stepwiselm(tbl, 'interactions', 'Upper', 'quadratic');
    mdlCoeff = [mdl.Coefficients.Estimate(1:end),mdl.Coefficients.SE(1:end),...
        mdl.Coefficients.tStat(1:end), mdl.Coefficients.pValue(1:end)];
    mdlR2Adj = mdl.Rsquared.Adjusted;
    alpha = 0.05;
    mdlCI = coefCI(mdl,alpha);                                              % lower and upper (1-alpha) confidence limits
    
    %% Validation
    DataVal = DataVal(xRangeVal,yRange);
    QnVal = polyval([mQ tQ], DataVal(:,1));
    pnVal = polyval([mp tp], DataVal(:,2));
    TnVal = polyval([mT tT], DataVal(:,3));
    Exp = DataVal(:, result+3);
    
    Pred = predict(mdl, [QnVal, pnVal, TnVal]);
    Error = (Pred - Exp)./Exp;
    
    %% plots
    if plotting == 1
        figure()
        plotDiagnostics(mdl,'cookd')
        
        figure()
        plotResiduals(mdl,'probability')
        
        figure()
        plotAdded(mdl)
        
        figure()
        h = bar(mdl.Coefficients.Estimate(2:end));
        set(h,'facecolor',[0.8 0.8 0.9])
        legend('Coefficient')
        set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
        set(gca,'xticklabel',mdl.CoefficientNames(2:end))
        ylabel(strcat(Response, {' '},Unit))
        xlabel('Normalized Coefficient')
        title('Quadratic Model Coefficients')
        
        figure()
        bar([Exp Pred])
    end
    
    %% create xlsx file
    filename = strcat('DoE_', modulename, '.xlsx');
    outname = fullfile(pwd, filename);
    warning('off','MATLAB:xlswrite:AddSheet');
    if exist(filename, 'file') == 0
        xlswrite(outname, [Variables;Units], 'DoE', 'A1')
        xlswrite(outname, Data, 'DoE', 'A3')
    end
    
    xlswrite(outname, cell(20,10),Response{1}, 'A1');
    xlswrite(outname, mdl.CoefficientNames', Response{1}, 'A2');
    xlswrite(outname, {'Estimate';'SE';'tStat';'pValue'; ...
        'CIlow';'CIhigh'}', Response{1}, 'B1');
    xlswrite(outname, [mdlCoeff, mdlCI], Response{1},'B2');
    xlswrite(outname, {'R2Adjusted'}, Response{1}, 'I1');
    xlswrite(outname, mdlR2Adj, Response{1}, 'I2');

    % delete empty sheets (Tabelle 1) in data file
    file = fullfile(pwd, filename);
    newExcel = actxserver('excel.application');
    excelWB = newExcel.Workbooks.Open(file,0,false);
    newExcel.Sheets.Item(1).Delete;
    newExcel.Sheets.Item(1).Delete;
    newExcel.Sheets.Item(1).Delete;
    excelWB.Save();
    excelWB.Close();
    newExcel.Quit();
    delete(newExcel);
    
    %% display model
    disp(mdl)
end