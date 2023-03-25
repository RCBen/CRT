if mix == 1
    HG = HGr;
    sample = sampler;
elseif mix == 2
    HG = HGm;
    sample = samplem;
end

S = repmat(sample,size(Theta,2),1);
Y = reshape(repmat(HG,1,size(Theta,2)),length(sample)*size(Theta,2),1);
if prop == 6 || prop == 9
    X = Theta(:) + 273.15;
    Z = reshape(eval(sheet{prop}),length(X),1);
else    
    X = reshape(repmat(Theta + 273.15,length(sample),1),length(sample)*size(Theta,2),1);
    Z = reshape(eval(sheet{prop}),length(sample)*size(Theta,2),1);
end

row = find(isnan(X));
S = S(~ismember(1:size(S,1), row));
X = X(~ismember(1:size(X,1), row));
Y = Y(~ismember(1:size(Y,1), row));
Z = Z(~ismember(1:size(Z,1), row));
    
if prop == 5 || prop == 6
    Z = log(Z);
    sheet{prop} = strcat('log',sheet{prop});
    if prop == 5
        if mix == 1
            model = 'poly33';
        elseif mix == 2
            model = 'poly32';
        end
    elseif prop == 6
        model = 'quadratic';
%         model = 'interactions';
    end
% elseif prop == 9
%     model = 'quadratic';
else
    model = 'interactions';
end

Expmt = table(S, ...
    X, Y, Z, 'VariableNames',{'sample','T','HG', sheet{prop}});

%% Saving to excel file
warning('off','MATLAB:xlswrite:AddSheet');

File = strcat(sheet{prop}, mixture{mix},'Fit.xlsx');
ExpmtSheet = strcat('Expmt.', method{hg});
mdlSheet = strcat('mdl.', sheet{prop}, '.', method{hg});

if ~exist(File, 'file') == 0
    [~, sheets] = xlsfinfo(File);
end
if exist(File, 'file') == 0 || ~any(strcmp(sheets, ExpmtSheet))
    writetable(Expmt, File, 'Sheet', ExpmtSheet)
    
    % delete empty sheets in data file
%     file = fullfile(pwd, File);
%     newExcel = actxserver('excel.application');
%     excelWB = newExcel.Workbooks.Open(file,0,false);
%     newExcel.Sheets.Item(1).Delete;
%     newExcel.Sheets.Item(1).Delete;
%     newExcel.Sheets.Item(1).Delete;
%     excelWB.Save();
%     excelWB.Close();
%     newExcel.Quit();
%     delete(newExcel);
else
    writecell(cell(30,1), File, 'Sheet', ExpmtSheet, 'Range', 'F1:F30')
    writecell(cell(19,7), File, 'Sheet', mdlSheet, 'Range', 'A2:G20')
%     xlswrite(File, cell(30,1), ExpmtSheet, 'F1:F30')
%     xlswrite(File, cell(19,7), mdlSheet, 'A2:G20')
end

%% Estimate the coefficients of the model using the fitlm 
% mdl2 = fitlm(Expmt(:,2:4), model, 'RobustOpts', 'on');
mdl2 = fitlm(Expmt(:,2:4),model);

%% Remove point with large Cook's distance from the model
larg = find((mdl2.Diagnostics.CooksDistance)>3*mean(mdl2.Diagnostics.CooksDistance));
Resultlarg = [S(larg), X(larg) - 273.15, Y(larg), Z(larg)];

if excl == 0
    larg = [];
    mdl = mdl2;
elseif excl == 1
    mdl = fitlm(Expmt(:,2:4),model, 'Exclude',larg);

%     if mdl2.Rsquared.Ordinary > 0.99 || mdl.Rsquared.Ordinary < mdl2.Rsquared.Ordinary + 0.05*(1-mdl2.Rsquared.Ordinary) 
%         larg = [];
%         mdl = mdl2;
%     else
    writecell({'excluded'}, File, 'Sheet', ExpmtSheet, 'Range', 'F1')
    writematrix(larg+1, File, 'Sheet', ExpmtSheet, 'Range', 'F2')
%     xlswrite(File, {'excluded'}, ExpmtSheet, 'F1')
%     xlswrite(File, larg+1, ExpmtSheet, 'F2')
%     end
end

mdlCoeff = [mdl.Coefficients.Estimate(1:end),mdl.Coefficients.SE(1:end),...
    mdl.Coefficients.tStat(1:end), mdl.Coefficients.pValue(1:end)];
mdlR2 = mdl.Rsquared.Ordinary;
alpha = 0.05;
mdlCI = coefCI(mdl,alpha);                                                  % lower and upper (1-alpha) confidence limits

% %% Diagnostics
% figure()
% plotDiagnostics(mdl,'cookd')                                                % Look for points with large Cook's distance.
% figure()
% plotResiduals(mdl,'probability')                                            % full-normal plot
% figure()
% plotAdded(mdl)                                                              % Create an added variable plot for the model as a whole.

%% Saving to excel file
writecell(mdl.CoefficientNames', File, 'Sheet', mdlSheet, 'Range', 'A2')
writecell({'Estimate';'SE';'tStat';'pValue'; 'CIlow';'CIhigh'}', File, 'Sheet', mdlSheet, 'Range', 'B1')
writematrix([mdlCoeff, mdlCI], File, 'Sheet', mdlSheet, 'Range', 'B2')
writecell({'Rsquared'}, File, 'Sheet', mdlSheet, 'Range', 'I1')
writematrix(mdlR2, File, 'Sheet', mdlSheet, 'Range', 'I2')
% xlswrite(File, mdl.CoefficientNames', mdlSheet, 'A2');
% xlswrite(File, {'Estimate';'SE';'tStat';'pValue'; 'CIlow';'CIhigh'}', mdlSheet, 'B1');
% xlswrite(File, [mdlCoeff, mdlCI], mdlSheet,'B2');
% xlswrite(File, {'Rsquared'}, mdlSheet, 'I1');
% xlswrite(File, mdlR2, mdlSheet, 'I2');

% %% Display the magnitudes of the coefficients (for normalized values) in a bar chart.
% figure()
% h = bar(mdl.Coefficients.Estimate(2:end));
% set(h,'facecolor',[0.8 0.8 0.9])
% legend('Coefficient')
% set(gcf,'units','normalized','position',[0.05 0.4 0.35 0.4])
% set(gca,'xticklabel',mdl.CoefficientNames(2:end))
% ylabel(strcat(sheet{prop}, ' (',unit{prop}, ')'))
% xlabel('Coefficient')
% title('Quadratic Model Coefficients')


%% plot model vs exp
mdlpred = predict(mdl, [X, Y]);

mdlprednew = NaN(size(mdlpred,1) + numel(row),1);                           % Preallocate
idx  = setdiff(1:size(mdlpred,1)+numel(row),row);                           % all positions of Anew except pos
mdlprednew(idx,:) = mdlpred;

figure()
surf(Theta + 273.15, HG, reshape(mdlprednew, length(sample),size(Theta,2)))
hold on
scatter3(X, Y, Z,'filled', 'MarkerFaceColor', 'k')
if isempty(larg) == 1
    scatter3(X((mdl2.Diagnostics.CooksDistance)>3*mean(mdl2.Diagnostics.CooksDistance)), ...
        Y((mdl2.Diagnostics.CooksDistance)>3*mean(mdl2.Diagnostics.CooksDistance)), ...
        Z((mdl2.Diagnostics.CooksDistance)>3*mean(mdl2.Diagnostics.CooksDistance)),'filled', 'MarkerFaceColor', 'g')
else
    scatter3(X(larg), Y(larg), Z(larg),'filled', 'MarkerFaceColor', 'r')
end
ylabel('HG (-)')
zlabel(strcat(sheet{prop}, ' (',unit{prop}, ')'))
xlabel('T (K)')

%% Andrade correlation (only T) / DIPPR 101
if prop == 5
    beta0 = ones(3,1)*10^-6;
    modelfun = 'logeta ~ A + B/T + C*log(T)';
    for i = 1:length(sample)        
        mdlA = fitnlm(Expmt(i:length(sample):end,2:4), modelfun, beta0);
        largA = find((mdlA.Diagnostics.CooksDistance)>3*mean(mdlA.Diagnostics.CooksDistance));
        
        if excl == 1
            writecell({'excluded'}, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'P1')
            writematrix(largA, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'P2')
%             xlswrite(File, {'excluded'}, strcat(mdlSheet,'.',num2str(sample(i))), 'P1')
%             xlswrite(File, largA, strcat(mdlSheet,'.',num2str(sample(i))), 'P2')
            mdlA = fitnlm(Expmt(i:length(sample):end,2:4), modelfun, beta0, 'Exclude', largA);
        end
        
        mdlACoeff = [mdlA.Coefficients.Estimate(1:end),mdlA.Coefficients.SE(1:end),...
                mdlA.Coefficients.tStat(1:end), mdlA.Coefficients.pValue(1:end)];
        mdlAR2 = mdlA.Rsquared.Ordinary;
        mdlACI = coefCI(mdlA,alpha);
        
        writecell(mdlA.CoefficientNames', File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'A2')
        writecell({'Estimate';'SE';'tStat';'pValue'; 'CIlow';'CIhigh'}', File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'B1')
        writematrix([mdlACoeff, mdlACI], File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'B2')
        writecell({'Rsquared'}, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'I1')
        writematrix(mdlAR2, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'I2')
%         xlswrite(File, mdlA.CoefficientNames', strcat(mdlSheet,'.',num2str(sample(i))), 'A2');
%         xlswrite(File, {'Estimate';'SE';'tStat';'pValue'; 'CIlow';'CIhigh'}', strcat(mdlSheet,'.',num2str(sample(i))), 'B1');
%         xlswrite(File, [mdlACoeff, mdlACI], strcat(mdlSheet,'.',num2str(sample(i))),'B2');
%         xlswrite(File, {'Rsquared'}, strcat(mdlSheet,'.',num2str(sample(i))), 'I1');
%         xlswrite(File, mdlAR2, strcat(mdlSheet,'.',num2str(sample(i))), 'I2');
    end  
elseif prop == 6
    beta0 = [100, -1000, -10];
    modelfun = 'logpLV ~ A + B/T + C*log(T)';
%     modelfun = 'logpLV ~ (A + B/T + C*log(T/298.15))/8.3143';               % modified Antoine
    options = optimoptions('fsolve','Display','none');
    ThetaLVall = zeros(length(sample),1);
    exitflagALL = zeros(length(sample),1);
    for i = 1:length(sample)
        [ThetaLVall(i), ~, exitflagALL(i)] = fsolve(@(T) ...
            mdl.Coefficients.Estimate(1)+ mdl.Coefficients.Estimate(2)*T + mdl.Coefficients.Estimate(3)*HG(i) + mdl.Coefficients.Estimate(4)*T*HG(i) + ...
            mdl.Coefficients.Estimate(5)*T.^2 +mdl.Coefficients.Estimate(6)*HG(i).^2 - log(1013.25), ...
            350+273.15, options);
        ThetaLVall(i) = ThetaLVall(i) - 273.15;
%         ThetaLVall(i) = fsolve(@(T) ...
%             mdl.Coefficients.Estimate(1)+ mdl.Coefficients.Estimate(2)*T + mdl.Coefficients.Estimate(3)*HG(i) + mdl.Coefficients.Estimate(4)*T*HG(i) - ...
%             log(1013.25), 350+273.15, options) - 273.15;
        
        XA = Theta(i,:) + 273.15;
        XA = XA(~ismember(1:size(XA,2), find(isnan(XA))));
        YA = repmat(HG(i,:),1,length(XA));
        ZA = log(pLV(i,:));
        ZA = ZA(~ismember(1:size(ZA,2), find(isnan(ZA))));

        ExpmtA = table(XA', YA', ZA', ....
            'VariableNames',{'T','HG', sheet{prop}});
            
        mdlA = fitnlm(ExpmtA, modelfun, beta0); 
        largA = find((mdlA.Diagnostics.CooksDistance)>3*mean(mdlA.Diagnostics.CooksDistance));
            
        if excl == 1
            writecell({'excluded'}, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'P1')
            writematrix(largA, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'P2')
%             xlswrite(File, {'excluded'}, strcat(mdlSheet,'.',num2str(sample(i))), 'P1')
%             xlswrite(File, largA, strcat(mdlSheet,'.',num2str(sample(i))), 'P2')
            mdlA = fitnlm(ExpmtA, modelfun, beta0, 'Exclude', largA);            
        end            

        mdlACoeff = [mdlA.Coefficients.Estimate(1:end),mdlA.Coefficients.SE(1:end),...
                mdlA.Coefficients.tStat(1:end), mdlA.Coefficients.pValue(1:end)];
        mdlAR2 = mdlA.Rsquared.Ordinary;
        mdlACI = coefCI(mdlA,alpha);
        
%         [mdlA.Coefficients.Estimate,~,residual] = lsqnonlin(@(p_est) p_est(1).*ones(length(XA'),1) + p_est(2)./XA' + p_est(3).*log(XA') -  ZA', ...
%                 [100, -1000, -10]', [0, -Inf, -Inf]', [Inf, Inf, Inf]', optimset( 'Algorithm', 'trust-region-reflective', 'Display', 'None'));
%             
%         mdlAR2 = 1-sum(residual.^2)./sum((ZA' - mean(ZA)).^2);
        
        writecell(mdlA.CoefficientNames', File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'A2')
        writecell({'Estimate';'SE';'tStat';'pValue'; 'CIlow';'CIhigh'}', File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'B1')
        writematrix([mdlACoeff, mdlACI], File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'B2')
        writecell({'Rsquared'}, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'I1')
        writematrix(mdlAR2, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'I2')
%         xlswrite(File, mdlA.CoefficientNames', strcat(mdlSheet,'.',num2str(sample(i))), 'A2');
%         xlswrite(File, {'Estimate';'SE';'tStat';'pValue'; 'CIlow';'CIhigh'}', strcat(mdlSheet,'.',num2str(sample(i))), 'B1');
%         xlswrite(File, [mdlACoeff, mdlACI], strcat(mdlSheet,'.',num2str(sample(i))),'B2');
%         xlswrite(File, {'Rsquared'}, strcat(mdlSheet,'.',num2str(sample(i))), 'I1');
%         xlswrite(File, mdlAR2, strcat(mdlSheet,'.',num2str(sample(i))), 'I2');    

        if i == 1
            ThetaLV = zeros(length(sample),1);
        end
        ThetaLV(i,1) = fsolve(@(T) ...
            mdlA.Coefficients.Estimate(1) + mdlA.Coefficients.Estimate(2)/T + mdlA.Coefficients.Estimate(3)*log(T) - log(1013.25), ...
            350+273.15, options) - 273.15;
            
        writecell({'TLV'}, File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'K1')
        writematrix(ThetaLV(i,1), File, 'Sheet', strcat(mdlSheet,'.',num2str(sample(i))), 'Range', 'K2')
%         xlswrite(File, {'TLV'}, strcat(mdlSheet,'.',num2str(sample(i))), 'K1');
%         xlswrite(File, ThetaLV(i,1), strcat(mdlSheet,'.',num2str(sample(i))), 'K2');
    end
    writecell({'TLV'}, File, 'Sheet', mdlSheet, 'Range', 'K1')
    writematrix(ThetaLVall, File, 'Sheet', mdlSheet, 'Range', 'K2')
    writecell({'exitflag'}, File, 'Sheet', mdlSheet, 'Range', 'L1')
    writematrix(exitflagALL, File, 'Sheet', mdlSheet, 'Range', 'L2')
%     xlswrite(File, {'TLV'}, mdlSheet, 'K1');
%     xlswrite(File, ThetaLVall, mdlSheet, 'K2');
end

%% save variables
% save(strcat(sheet{prop},'Fit'))
    
% end