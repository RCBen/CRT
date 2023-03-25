clear variables
if isempty(which('getConstants'))
    mydir = pwd;
    idcs = strfind(mydir,'\');
    p.subdir = mydir(1:idcs(end-2)-1);
    addpath(strcat(p.subdir,'\functions'))
end
const = getConstants;

%% Definitions                                                                 
p.smoothRaw = 1;                                                            % smooth raw data    0 off
                                                                            %                    1 on
p.smoothResults = 0;                                                        % smooth results     0 off
                                                                            %                    1 on
p.DPt = 0;                                                                  % dispersion for mPt 1 on
                                                                            %                    0 off                                                                           
resultfileX = 'Ergebnisse.xlsx';
resultfile ='EvaluationParticle.mat';
smoothmethod = 'rlowess';                                                                            

% reactor                                                              
pr = 1;                                                                     % reaction pressure (bar)   
nP = numel(pr);

Thetarstr = {'310' '300' '290' '280'};                            
Thetar = str2double(Thetarstr);                                             % reaction temperature (°C)
nTheta = numel(Thetar);
Tr = Thetar + 273.15;                                                       % reaction temperatures (K)

Vr = 100;                                                                   % reactor volume (ml)

% LOHC
HG0 = 0.95;                                                                 % initial hydrogenation grade (-)
mL0 = 29.1;                                                                 % initial amount of liquid phase (g);
nL0 = mL0/(HG0*const.M18 + (1-HG0)*const.M0);                               % initial molar amount of liquid phase (mol);

nH2max = mL0*const.wH218/const.MH2;                                         % maximum of releasable H2 regarding to mL0(i. e. HG = 1) (mol)
nH20 = HG0*nH2max;                                                          % maximum of releasable H2 regarding to HG0 (mol)

% catalyst
r = const.r;                                                                % radius of uncoated thermocouple (m)

dstr = {'20' '40' '60' '80' '110' '150' '200'};
% dRange = [1:5,7];
dRange = 1:7;
dstr = dstr(dRange);
dH = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'd', ...                % thickness of active layer for thermocouple (m)
    'Range', 'C17:C23')*1e-06;            
dH = dH(dRange);
nd = numel(dH);
l = 0.020;                                                                  % length of active layer for thermocouple (m)
 
mPt = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'mPt', ...             % mass of Pt for thermocouples & thermocouples and spheres(g)  
    'Range', 'B12:C18')/1000; 
if p.DPt == 1
    DPt = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'CO', ...         % dispersion of Pt for thermocouples & spheres (-)  
        'Range', 'B27:C33');
    DPt(DPt>1)=1;
    mPt = mPt.*DPt;
end
mPt = mPt(dRange,:);
H.xPt = (mPt(:,1)/const.MPt)/nL0;                                           % molar fraction of Pt regarding to nL0 for thermocouples (-)
S.xPt = (sum(mPt,2)/const.MPt)/nL0;                                         % molar fraction of Pt regarding to nL0 for thermocouples and spheres (-)

wPtLH = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'mPt', ...          % mass fraction of Pt regarding to mass of active layer for thermocouples (-)
    'Range', 'C2:C8'); 
wPtLH = wPtLH(dRange,1);
wPtLS = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'mPt', ...          % mass fraction of Pt regarding to mass of active layer for spheres (-)
    'Range', 'F2:F8'); 
wPtLS = wPtLS(dRange,1);

dS = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'd', ...                % thickness of active layer for spheres (m)
    'Range', 'C4:C10')*1e-06;   
dS = dS(dRange);
mc = readmatrix(fullfile(pwd, resultfileX), 'Sheet', 'mPt', ...              % mass of spheres (g)
    'Range', 'D2:D8');  
mc = mc(dRange,1);
H.mc = (mPt(:,1)./mPt(:,2)).*mc;                                            % mass of cat for thermocouples (g)
S.mc = H.mc + mc;                                                           % mass of cat for thermocouples and spheres(g)

% S.dpore80 = 3.466*1e-09;                                                    % average pore diameter for 80 um layer of spheres (m)
% H.dpore80 = 13.73*1e-09;                                                    % average pore diameter for 80 um layer of thermocouple (m)

dHr = 1000*getDHr(const, pr, Tr)/9;                                         % reaction enthalpy for dehydrogenation of H18-DBT (J molH2^-1)

%% Read in raw data (T, FTC)
startTime = tic;

Cat = {'H' 'S'};                                                            % type of catalyst (-)
Variables = {'t' 'Tcat' 'TLOHC' 'xH2'};
nVariables = numel(Variables);
yRange = [1, 3:5];
yt = yRange(sym(Variables) == 't');
yTcat = yRange(sym(Variables) == 'Tcat');
yxH2 = yRange(sym(Variables) == 'xH2');

pathname = fullfile(pwd, 'Auswertung');

Thetamin = 250;                                                             % min temperature for dehydrogenation of H18-DBT (°C)
xRangeHmin = ones(nd,1);
xRangeSmin = ones(nd,nTheta);
for i = 1:nd
    % thermocouples
    filename = strcat(Cat{1}, dstr{i}, '.xlsx');
    if exist(fullfile(pathname,filename), 'file')
        for j = 1
            Data = readmatrix(fullfile(pathname, filename), 'Sheet', ...
                Thetarstr{j}, 'Range', 'A2');            
            if ~isempty(Data(~isnan(Data))) && ~all(isnan(Data(:,yTcat)))
                xRangeHmin(i,j) = find(Data(:,yTcat) > Thetamin, 1, 'first')-1;
                Data(:, yxH2) = Data(:, yxH2)/100;                          % convert xH2 from [%] to [-]
                H.Baseline(i,j) = ...
                    mean(Data(~isoutlier(Data(1:xRangeHmin(i,j),yxH2)),yxH2));
                Data(:, yxH2) = Data(:, yxH2) - H.Baseline(i,j);            % delete baseline
                Data(Data(:, yxH2)<0, yxH2) = 0;
                
                Data = Data(xRangeHmin(i,j):end,:);                         % erase heating phase until Tmin
                Data(:, yt) = Data(:, yt) - Data(1, yt);                    % set t(Theat = Tmin) = 0                
                for k = 1:nVariables
                    H.(Variables{k}){i,j} = Data(:, yRange(k));
                end
            end
        end
    end

    % thermocouples + spheres
    filename = strcat(Cat{2}, dstr{i}, '.xlsx');
    if exist(fullfile(pathname,filename), 'file')
        for j = 1:nTheta
            Data = readmatrix(fullfile(pathname, filename), 'Sheet', ...
                Thetarstr{j}, 'Range', 'A2'); 
            if ~isempty(Data(~isnan(Data))) && ~all(isnan(Data(:,yTcat)))
                xRangeSmin(i,j) = find(Data(:,yTcat) > Thetamin, 1, 'first')-1;
                Data(:, yxH2) = Data(:, yxH2)/100;                          % convert xH2 from [%] to [-]
                S.Baseline(i,j) = ...
                    mean(Data(~isoutlier(Data(1:xRangeSmin(i,j),yxH2)),yxH2));
                Data(:, yxH2) = Data(:, yxH2) - S.Baseline(i,j);            % delete baseline
                Data(Data(:, yxH2)<0, yxH2) = 0;
                
                Data = Data(xRangeSmin(i,j):end,:);                         % erase heating phase until Tmin
                Data(:, yt) = Data(:, yt) - Data(1, yt);                    % set t(Theat = Tmin) = 0                
                for k = 1:nVariables
                    S.(Variables{k}){i,j} = Data(:, yRange(k));
                end
            end
        end
    end
end
calcTimeRead = toc(startTime);

%% Process data
QN2 = 0.5;                                                                  % volume flow rate of carrier gas NS (l min^-1)
M = 30;                                                                     % start entry for cumtrapz
MT = 721;                                                                   % number of entries for Q used to determine max(Q) 
Qmaxvalue = 1000;
for i = 1:nd
    % thermocouples
    for j = 1
        if ~isempty(H.t{i,j})
            H.dTRaw{i,j} = H.TLOHC{i,j} - H.Tcat{i,j};
            if p.smoothRaw == 1
                H.dT{i,j} = smooth(H.t{i,j}, H.dTRaw{i,j}, smoothmethod);
            else
                H.dT{i,j} = H.dTRaw{i,j};
            end

            H.QH2tRaw{i,j} = H.xH2{i,j}*QN2./(1 - H.xH2{i,j});              % temporal QH2 at t (l min^-1)
            if p.smoothRaw == 1
                H.QH2t{i,j} = smooth(H.t{i,j}, H.QH2tRaw{i,j}, smoothmethod);
            else 
                H.QH2t{i,j} = H.QH2tRaw{i,j};
            end
            H.QH2m{i,j} = cumtrapz(H.t{i,j}(M:end), H.QH2t{i,j}(M:end))./...
                (H.t{i,j}(M:end)-H.t{i,j}(M));
            H.QH2m{i,j} = [H.QH2t{i,j}(1:M-1); H.QH2m{i,j}];                % mean QH2 from t=0 to t (l min^-1)
            
            H.Prod{i,j} = H.QH2t{i,j}*const.rhoH2n/mPt(i,1);                % temporal productivity at t (gH2 gPt^-1 min^-1)
            H.Prodm{i,j} = H.QH2m{i,j}*const.rhoH2n/mPt(i,1);               % mean productivity from t=0 to t (gH2 gPt^-1 min^-1)
            H.Pr{i,j} = H.QH2t{i,j}/60*const.rhoH2n./const.MH2*dHr(j);
            H.LambdaEff{i,j} = H.Pr{i,j}./ ...
                (2*pi*l/log((const.r+dH(i))/const.r)*H.dT{i,j});
            
            H.HG{i,j} = (nH20 - ...
                (cumtrapz(H.t{i,j}, H.QH2t{i,j})*const.rhoH2n/const.MH2))./ ...
                nH2max;
            if ~isempty(H.HG{i,j})
                H.tend(i,j) = H.t{i,j}(end);
                H.HGend(i,j) = H.HG{i,j}(end);                
            end
            H.X{i,j} = 1-H.HG{i,j}./HG0;
            H.Q{i,j} = [H.X{i,j}(1); ...                                    % local heat flow from reaction at t (W) 
                H.X{i,j}(2:end-1)-H.X{i,j}(1:end-2); H.X{i,j}(end)]* ...
                dHr(j)*9*nL0./(H.t{i,j}(2)-H.t{i,j}(1));
            H.Qm{i,j} = H.X{i,j}*dHr(j)*9*nL0./(H.t{i,j}-H.t{i,j}(1));      % mean heat flow from reaction from t=0 to t (W)
            if size(H.t{i,j},1) >= MT
                H.Qmax{i,j} = max(H.Q{i,j}(H.Q{i,j}(1:MT) < Qmaxvalue));
                H.Qmmax{i,j} = max(H.Qm{i,j}(H.Qm{i,j}(1:MT) < Qmaxvalue));
            else
                H.Qmax{i,j} =  max(H.Q{i,j}(H.Q{i,j} < Qmaxvalue));
                H.Qmmax{i,j} =  max(H.Q{i,j}(H.Qm{i,j} < Qmaxvalue));
            end
        end
    end

    % thermocouples + spheres
    for j = 1:nTheta
        if ~isempty(S.t{i,j})
            S.dTRaw{i,j} = S.TLOHC{i,j} - S.Tcat{i,j};
            if p.smoothRaw == 1 && ~(i == 3 && j == 1)
                S.dT{i,j} = smooth(S.t{i,j}, S.dTRaw{i,j}, smoothmethod);
            else
                S.dT{i,j} = S.dTRaw{i,j};
            end
            
            S.QH2tRaw{i,j} = S.xH2{i,j}*QN2./(1 - S.xH2{i,j});              % temporal QH2 at t (l min^-1)
            if p.smoothRaw == 1 && ~(i == 3 && j == 1)
                S.QH2t{i,j} = smooth(S.t{i,j}, S.QH2tRaw{i,j}, smoothmethod);
            else 
                S.QH2t{i,j} = S.QH2tRaw{i,j};    
            end
            S.QH2m{i,j} = cumtrapz(S.t{i,j}(M:end), S.QH2t{i,j}(M:end))./... 
                (S.t{i,j}(M:end)-S.t{i,j}(M));
            S.QH2m{i,j} = [S.QH2t{i,j}(1:M-1); S.QH2m{i,j}];                % mean QH2 from t=0 to t (l min^-1)
            
            S.Prod{i,j} = S.QH2t{i,j}*const.rhoH2n/sum(mPt(i,:),2);         % temporal productivity at t (gH2 gPt^-1 min^-1)
            S.Prodm{i,j} = S.QH2m{i,j}*const.rhoH2n/sum(mPt(i,:),2);        % mean productivity from t=0 to t (gH2 gPt^-1 min^-1)
            S.Pr{i,j} = S.QH2t{i,j}/60*const.rhoH2n./const.MH2*dHr(j);          
            S.LambdaEff{i,j} = S.Pr{i,j}*(mPt(i,1)/sum(mPt(i,:),2))./ ...
                (2*pi*l/log((const.r+dS(i))/const.r)*S.dT{i,j});

            S.HG{i,j} = (nH20 - ...
                (cumtrapz(S.t{i,j}, S.QH2t{i,j})*const.rhoH2n/const.MH2))./ ...
                nH2max;
            if ~isempty(S.HG{i,j})
                S.tend(i,j) = S.t{i,j}(end);
                S.HGend(i,j) = S.HG{i,j}(end);
            end
            S.X{i,j} = 1-S.HG{i,j}./HG0;
            S.Q{i,j} = [S.X{i,j}(1); ...                                    % local heat flow from reaction at t (W) 
                S.X{i,j}(2:end-1)-S.X{i,j}(1:end-2); S.X{i,j}(end)]* ...
                dHr(j)*9*nL0./(S.t{i,j}(2)-S.t{i,j}(1));
            S.Qm{i,j} = S.X{i,j}*dHr(j)*9*nL0./(S.t{i,j}-S.t{i,j}(1));      % mean heat flow from reaction from t=0 to t (W)   
            if size(S.t{i,j},1) >= MT
                S.Qmax{i,j} = max(S.Q{i,j}(S.Q{i,j}(1:MT) < Qmaxvalue));
                S.Qmmax{i,j} = max(S.Qm{i,j}(S.Qm{i,j}(1:MT) < Qmaxvalue));
            else
                S.Qmax{i,j} =  max(S.Q{i,j}(S.Q{i,j} < Qmaxvalue));
                S.Qmmax{i,j} =  max(S.Q{i,j}(S.Qm{i,j} < Qmaxvalue));
            end
        end
    end
end
calcTimeProcess = toc(startTime) - calcTimeRead;

%% mass balance
p.pr = pr;                                                                  % reaction pressure (bar)  
p.Vr = Vr;                                                                  % reactor volume (ml)
p.nL0 = nL0;                                                                % initial molar amount of liquid phase (mol)

for i = 1:nd
    % thermocouples
    for j = 1   
        if ~isempty(H.t{i,j}) && issorted(H.HG{i,j}, 'monotonic')
            p.Tr = Tr(j);                                                   % reaction temperature (K)
            p.mc = H.mc(i);                                                 % mass of cat (g)
            p.t = H.t{i,j}*60;                                              % reaction time (s)
            p.HG = H.HG{i,j};                                               % hydrogenation grade (-)
            H.mb{i,j} = myMassBalance(p);
        end
    end
    
    % thermocouples + spheres
    for j = 1:nTheta
        if ~isempty(S.t{i,j}) && issorted(S.HG{i,j}, 'monotonic')
            p.Tr = Tr(j);                                                   % reaction temperature (K)
            p.mc = S.mc(i);                                                 % mass of cat (g)
            p.t = S.t{i,j}*60;                                              % reaction time (s)
            p.HG = S.HG{i,j};                                               % hydrogenation grade (-)
            S.mb{i,j} = myMassBalance(p);
        end
    end
end
p = rmfield(p, 'HG');
calcTimeMB = toc(startTime) - calcTimeProcess;

%% Results @HG
Var = [Variables, 'dT' 'QH2m' 'Prod' 'Pr' 'LambdaEff' 'HG' ];
nVar = numel(Var);
Units = {'min' '°C' '°C' '-' '°C' 'Nml min^{-1}' 'gH2 gPt^{-1} min^{-1}' ...
    'W' 'W m^{-1} K^{-1}' '-'};

dHG = 0.05;
HGend = 0.60;
HG = (HG0-dHG):-dHG:HGend;
nHG = numel(HG);

span = 7;
N = (span-1)/2;                                                               % neighbouring values                                                                         

H.HGend(H.HGend < HGend & H.HGend ~= 0) = HGend;
S.HGend(S.HGend < HGend & S.HGend ~= 0) = HGend;
H.idxHG = cell(nd,nTheta);
S.idxHG = cell(nd,nTheta);
H.idxHGrange = cell(nd,nTheta,nHG);
S.idxHGrange = cell(nd,nTheta,nHG);
for i = 1:nd
    % thermocouples
    for j = 1        
        if ~isempty(H.t{i,j})
            for h = 1:nHG
                helpH = H.HG{i,j}(~isnan(H.HG{i,j}));
                if helpH(end) <= HG(h)
                    H.idxHG{i,j}(1,h) = find(H.HG{i,j} >= HG(h), 1, 'last');
                    H.idxHGrange{i,j,h} = ((H.idxHG{i,j}(1,h)-N):1:(H.idxHG{i,j}(1,h)+N))';
                end
            end
            
            H.idxProd{i,j} = find(H.Prod{i,j} == max(H.Prod{i,j}), 1, 'first');
            H.idxProdrange{i,j} = ((H.idxProd{i,j}-N):1:(H.idxProd{i,j}+N))';
            H.idxProdrange{i,j} = H.idxProdrange{i,j}(H.idxProdrange{i,j}>0 & ...
                H.idxProdrange{i,j}<=numel(H.Prod{i,j}));
            
            Results.H.t{i,j} = H.t{i,j}([H.idxProd{i,j} H.idxHG{i,j}]);
            for k = (nVariables+1):nVar
                if p.smoothResults == 1 && ~(i == 3 && j == 1) && k < nVar
                    H.(Var{k}){i,j} = ...
                        smooth(H.t{i,j}, H.(Var{k}){i,j}, smoothmethod);
                end

                Res = H.(Var{k}){i,j}(H.idxProdrange{i,j});
                Results.H.(Var{k}){i,j}(1,1) = ...
                     mean(Res(~isnan(Res) & ~isoutlier(Res) & ~isinf(Res)));
                for h = 1:nHG
                    if ~isempty(H.idxHGrange{i,j,h})
                        Res = H.(Var{k}){i,j}(H.idxHGrange{i,j,h});
                        Results.H.(Var{k}){i,j}(h+1,1) = ...
                            mean(Res(~isnan(Res) & ~isoutlier(Res) & ~isinf(Res)));
                    end                    
                end
            end
        end
    end
    
    % thermocouples + spheres
    for j = 1:nTheta        
        if ~isempty(S.t{i,j})
            for h = 1:nHG
                helpS = S.HG{i,j}(~isnan(S.HG{i,j}));
                if helpS(end) <= HG(h)
                    S.idxHG{i,j}(1,h) = find(S.HG{i,j} >= HG(h), 1, 'last');
                     S.idxHGrange{i,j,h} = ((S.idxHG{i,j}(1,h)-N):1:(S.idxHG{i,j}(1,h)+N))';
                end
            end
            
            S.idxProd{i,j} = find(S.Prod{i,j} == max(S.Prod{i,j}), 1, 'first');
            S.idxProdrange{i,j} = ((S.idxProd{i,j}-N):1:(S.idxProd{i,j}+N))';
            S.idxProdrange{i,j} = S.idxProdrange{i,j}(S.idxProdrange{i,j}>0 & ...
                S.idxProdrange{i,j}<=numel(S.Prod{i,j}));

            Results.S.t{i,j} = S.t{i,j}([S.idxProd{i,j} S.idxHG{i,j}]);
            for k = (nVariables+1):nVar
                if p.smoothResults == 1 && ~(i == 3 && j == 1) && k < nVar
                    S.(Var{k}){i,j} = ...
                        smooth(S.t{i,j}, S.(Var{k}){i,j}, smoothmethod);
                end

                Res = S.(Var{k}){i,j}(S.idxProdrange{i,j});
                Results.S.(Var{k}){i,j}(1,1) = ...
                     mean(Res(~isnan(Res) & ~isoutlier(Res) & ~isinf(Res)));
                for h = 1:nHG
                    if ~isempty(S.idxHGrange{i,j,h})
                        Res = S.(Var{k}){i,j}(S.idxHGrange{i,j,h});
                        Results.S.(Var{k}){i,j}(h+1,1) = ...
                            mean(Res(~isnan(Res) & ~isoutlier(Res) & ~isinf(Res)));
                    end                    
                end
            end
            if ~isempty(S.t{i,j}) && ~isempty(S.idxHGrange{i,j, 1}) && numel(S.t{i,j}) >= S.idxHGrange{i,j,1}(end)
                Start = S.idxHG{i,j}(abs(HG - 0.9)<1e-6);
                helpSt = S.t{i,j}(~isnan(S.t{i,j}));
                helpSHG = S.HG{i,j}(~isnan(S.HG{i,j}));
                helpSQH2t = S.QH2t{i,j}(~isnan(S.QH2t{i,j}));
                if numel(S.idxHG{i,j})>=5 && numel(S.t{i,j}) >= S.idxHG{i,j}(5)
                    Ende = S.idxHG{i,j}(abs(HG - 0.7)<1e-6);
                else
%                     Ende = numel(S.t{i,j});
                    Ende =  numel(helpSt);
                end
%                 Results.S.Prodm9070(i,j) = trapz(S.t{i,j}( Start:Ende), S.QH2t{i,j}( Start:Ende))./... 
%                         (S.t{i,j}(Ende)-S.t{i,j}( Start))*const.rhoH2n/sum(mPt(i,:),2);
%                 Results.S.HG9070(i,j) = S.HG{i,j}(Ende);
                Results.S.Prodm9070(i,j) = trapz(helpSt( Start:Ende), helpSQH2t( Start:Ende))./...
                        (helpSt(Ende)-helpSt( Start))*const.rhoH2n/sum(mPt(i,:),2);           
                Results.S.HG9070(i,j) = helpSHG(Ende);
    %             Results.S.Prodm9070(i,j) = trapz(S.t{i,j}( S.idxHG{i,j}(abs(HG - 0.9)<1e-6):S.idxHG{i,j}(abs(HG - 0.7)<1e-6)), S.QH2t{i,j}( S.idxHG{i,j}(abs(HG - 0.9)<1e-6):S.idxHG{i,j}(abs(HG - 0.7)<1e-6)))./... 
    %                 (S.t{i,j}(S.idxHG{i,j}(abs(HG - 0.7)<1e-6))-S.t{i,j}( S.idxHG{i,j}(abs(HG - 0.9)<1e-6)))*const.rhoH2n/sum(mPt(i,:),2);
            end
        end
    end
end
calcTimeResult = toc(startTime) - calcTimeMB;

%% create xlsx file
warning('off','MATLAB:xlswrite:AddSheet');
h = struct2cell(Results.H);
s = struct2cell(Results.S);
for i = 1:nd
    % thermocouples
    for j = 1
        sheetname = strcat(Cat{1}, dstr{i});
        if ~isempty(Results.H.t{i,j})
            range = ['A', num2str((j-1)*11+3)];
            writematrix([Results.H.HG{i,j} cell2mat(arrayfun(@(x) ...
                h{x}{i,j}, 1:numel(h)-1, 'UniformOutput', false))], ...
                fullfile(pwd, resultfileX), 'Sheet', sheetname, ...
                'Range', range);
        end
    end
    
    % thermocouples + spheres
    for j = 1:nTheta
        sheetname = strcat(Cat{2}, dstr{i});       
        if ~isempty(Results.S.t{i,j})
            range = ['A', num2str((j-1)*11+3)];
            writematrix([Results.S.HG{i,j} cell2mat(arrayfun(@(x) ...
                s{x}{i,j}, 1:numel(s)-1, 'UniformOutput', false))], ...
                fullfile(pwd, resultfileX), 'Sheet', sheetname, ...
                'Range', range);
        end
    end
end
calcTimeWrite = toc(startTime) - calcTimeResult;

%% save variables
save(resultfile)