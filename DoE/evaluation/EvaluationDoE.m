clear variables;
% if ~exist('p', 'var') || ~isfield(p,'subdir')
    mydir = pwd;
    idcs = strfind(mydir,'\');
    p.subdir = mydir(1:idcs(end-2)-1);
    addpath(strcat(p.subdir,'\functions'))
% end
const = getConstants;
dt = 0.5;                                                                   % dt in min

%% Auswählen des Hauptordners - hier müssen alle Ergebnisdateien liegen
mainFolder = uigetdir;
[~,FolderName,~]=fileparts(mainFolder);
File = dir(fullfile(mainFolder,'*.txt'));
QLOHCsoll = str2double(File.name((3:4)+length(FolderName)));
Tsoll = str2double(File.name((7:9)+length(FolderName)));
p1soll = str2double(File.name((12:14)+length(FolderName)));
p2soll = str2double(File.name((17:19)+length(FolderName)));

GUI;
uiwait(msgbox('You need to use the point as decimal mark ','Note'));
waitfor(GUI);

%% Einlesen der Dateien
file = dir(fullfile(mainFolder, '*.dat'));
A = importfileDAT(fullfile (mainFolder, file.name));

file = dir(fullfile(mainFolder, '*.txt'));
T = importfileTXT(fullfile (mainFolder, file.name),3,inf);                  

file = dir(fullfile(mainFolder, '*.csv'));
P = importfileCSV(fullfile (mainFolder, file.name),1,inf);                  

if abs(size(P,1)-size(T,1)) > abs(size(P,1)-size(A,1))
    P = P([1,2:60*dt/A(3,1):end],:);
else
    P = P([1,2:60*dt/T(2,1):end],:);
end
A = A([1,2:60*dt/A(3,1):end],:);
T = T([1,2:60*dt/T(2,1):end],:);
P = P(3:end,:);                                                             % 3 korrigiert den verzögerten Start

% gleiche Länge
ix = min([length(A) length(T) length(P)]);
A = A(1:ix,:);
P = P(1:ix,:);
T = T(1:ix,:);

A(:,1) = A(:,1)/60;                                                         % Umrechnung in Minuten
T(:,1) = A(:,1);

%% Erzeugung der Matrizen
% D: Matrix mit allen Temperaturen
% E: Matrix mit allen Leistungen
% F: Matrix mit Pumpenzahlen, Drücken und Volumenströmen
D = [A(:,1) A(:,7) T(:,2:31) A(:,5) A(:,9) A(:,14) A(:,11) A(:,16) ...
    T(:,32:41) A(:,6) A(:,10) A(:,15) A(:,12) A(:,17)];
DHeader = { 'Zeit [min]' 'TIC-0 [°C]' 'A1 [°C]' 'A2 [°C]' 'A4 [°C]' ...
    'A6 [°C]' 'A7 [°C]' 'A8 [°C]' 'A9 [°C]' 'A10 [°C]' 'A12 [°C]' 'A14 [°C]'...
    'B6 [°C]' 'B7 [°C]' 'B8 [°C]' 'B9 [°C]' 'B10 [°C]' 'C6 [°C]' ...
    'C8 [°C]' 'C9 [°C]' 'C10 [°C]' 'D2 [°C]' 'D4 [°C]' 'D7 [°C]'...
    'D9 [°C]' 'D12 [°C]' 'D14 [°C]' 'E1 [°C]' 'E2 [°C]' 'E3 [°C]' ...
    'E4 [°C]' 'E5 [°C]' 'E6 [°C]' 'E7 [°C]' 'E8 [°C]' 'E9 [°C]'...
    'E10 [°C]' 'E11 [°C]' 'E12 [°C]' 'E13 [°C]' 'E14 [°C]' 'E15 [°C]' ...
    'F6 [°C]' 'F7 [°C]' 'F8 [°C]' 'F9 [°C]' 'F10 [°C]'... 
    'G6 [°C]' 'G7 [°C]' 'G8 [°C]' 'G9 [°C]' 'G10 [°C]' };

E = [A(:,1) P(:,1:6)];
EHeader = { 'Zeit [min]' 'PTIC-6 [W]' 'PTIC-7 [W]' 'PTIC-8 [W]' ...
    'PTIC-9 [W]' 'PTIC-10 [W]' 'PTIC-0 [W]'};
if ~isempty(find(E(:,2:6)>650, 1))                                          % suchen aller Leistungswerte > 650W
    disp('error')
end

F = [A(:,1:3) A(:,8) A(:,13) A(:,18)];
FHeader = { 'Zeit [min]' 'n [rpm]' 'MFM-H2 [l/min]' 'PIR-1 [bar]' ...
    'PIR-2 [bar]' 'AIR-CO [ppmv]' };

%% Erstellen einer Excel-Datei
[filename, pathname] = uiputfile(fullfile(mainFolder,strcat(FolderName,'.xlsx')), 'Choose a file name');
outname = fullfile(pathname, filename);
warning('off','MATLAB:xlswrite:AddSheet');

xlswrite(outname, D, 'T');
xlswrite(outname, DHeader, 'T');

xlswrite(outname, E, 'P');
xlswrite(outname, EHeader, 'P');

xlswrite(outname, F, 'n,VH2,p,CO');
xlswrite(outname, FHeader, 'n,VH2,p,CO');

% leere Reiter in Excel-Datei löschen
file = fullfile(pathname, filename);
newExcel = actxserver('excel.application');
excelWB = newExcel.Workbooks.Open(file,0,false);
newExcel.Sheets.Item(1).Delete;
newExcel.Sheets.Item(1).Delete;
newExcel.Sheets.Item(1).Delete;
excelWB.Save();
excelWB.Close();
newExcel.Quit();
delete(newExcel);

%% Berechnungen
% Headerline löschen
A = A(2:end,:);
P = P(2:end,:);
T = T(2:end,:);

D = D(2:end,:);
E = E(2:end,:);
F = F(2:end,:);

% Bestimmen der Ausgabezeiten
Tstart = 65;
[row, col] = find(D(:,33:37) >= Tstart);
tdTA = min(row);

tdTAoutput =  D(tdTA,1);
tdTEoutput = tdTE;
tSAoutput = tSA;
tSEoutput = tSE;

% Anpassen der Zeit auf dt-Intervall
if mod(tdTE,dt) < dt/2
    tdTE = (tdTE-mod(tdTE,dt))/dt + 1;
else 
    tdTE = (tdTE-mod(tdTE,dt))/dt + 2;
end

if mod(tSA,dt) < dt/2
    tSA = (tSA-mod(tSA,dt))/dt + 1;
else 
    tSA = (tSA-mod(tSA,dt))/dt + 2;
end

if mod(tSE,dt) < dt/2
    tSE = (tSE-mod(tSE,dt))/dt + 1;
else 
    tSE = (tSE-mod(tSE,dt))/dt + 2;
end

% DoE-Parameter
QLOHC = (mean(F(tSA:tSE,2)) + 7.287)/84.582;
TE = mean2(D(tSA:tSE,28:42));
sigmaTE = std2(D(tSA:tSE,28:42));
% TDE = mean2(D(tSA:tSE,22:42));
% sigmaTDE = std2(D(tSA:tSE,22:42));
p1 = mean(F(tSA:tSE,5));
p2 = 0;

% Wasserstoffvolumenstrom (l min^-1)
[~, cpH2] = getCp(const, [], [const.Tn, (TAbgas + 273.15)]+50, 1);
C = cpH2(1)*(TAbgas + 273.15)/(cpH2(2)*const.Tn);
QH2 = mean(F(tSA:tSE,3))*C;
sigmaQH2 = std(F(tSA:tSE,3));

% Produktivität (gH2 gPt^-1 min^-1)
Prod = QH2*const.rhoH2n/(mPt/1000);

% thermische Leistung Wasserstoff (W)
PthH2 = QH2/60*const.rhoH2n*const.LHVH2;                                                   

% idealer Leistungsbedarf (W)
% Reaktion
dHr = 1000*getDHr(const, p1, TE + 273.15)/9;                                % Reaktionsenthalpie (J molH2^-1)
Pr = QH2/60*const.rhoH2n/(2*1.00794)*dHr;

% Vorheizung
Tr = (const.Tn + (TE + 273.15))/2;
dT = TE + 273.15 - const.Tn;

rho = getRho(const, [], Tr, 1);
cp = getCp(const, [], Tr, 1)*1000/const.M18;
P0 = QLOHC/(60*1000)*rho*cp*dT;

% realer Leistungsbedarf (W)
% ohne Vorheizung
Pth = mean(sum(E(tSA:tSE,2:6),2));

% mit Vorheizung
Pthw0 = mean(sum(E(tSA:tSE,2:7),2));

% Wirkungsgrade (-)
% ohne Vorheizung
eta = PthH2/Pth;
etaId = PthH2/Pr;
Eta = eta/etaId;

% mit Vorheizung
etaw0 = PthH2/Pthw0;
etaIdw0 = PthH2/(Pr + P0);
Etaw0 = etaw0/etaIdw0;

% Energiebedarf Aufheizen (kJ)
t = A(tdTA:tdTE,1);
dtdT = A(tdTE) - A(tdTA);
EdT = trapz(60*t, sum(E(tdTA:tdTE,2:6),2))/1000;

%% Ergebnisse in Excel schreiben
xcoord = [1 1 1 2 2 2 2 2 3 3 2 2 2 2 2 2 2 2 2 1 1 2 2 3 3 1 1 1 1 1 ...
    2 2 2 2 2 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2]';
ycoord = [5 4 2 5 4 3 2 1 4 2 5 4 3 2 1 5 3 2 1 4 2 4 2 4 2 5 4 3 2 1 ...
    5 4 3 2 1 5 4 3 2 1 5 4 3 2 1 5 4 3 2 1]';
zcoord = [7 7 7 7 7 7 7 7 7 7 6 6 6 6 6 5 5 5 5 4 4 4 4 4 4 3 3 3 3 3 ...
    3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1]';

xlswrite(outname, {'x' 'y' 'z' 'T [°C]' 'Bezeichnung'}, 'T gemittelt','A1');
xlswrite(outname, xcoord, 'T gemittelt','A2');
xlswrite(outname, ycoord, 'T gemittelt','B2');
xlswrite(outname, zcoord, 'T gemittelt','C2');
xlswrite(outname, mean(D(tSA:tSE,3:end),1)', 'T gemittelt','D2');
xlswrite(outname, DHeader(3:end)', 'T gemittelt','E2');

Para1Header = {'QLOHC' 'p1' 'T' 'sigmaT'; '[ml min^-1]' '[bara]' '[°C]' '[°C]'};
Para1 = [QLOHC p1 TE sigmaTE];

Para2Header = {'TAbgas' 'C' 'QH2' 'sigmaQH2' 'Prod'; ...
    '[°C]' '[-]' '[Nml min^-1]' '[Nml min^-1]' '[gH2 gPt^-1 min^-1]'};
Para2 = [TAbgas C 1000*QH2 1000*sigmaQH2 Prod];

Para3Header = {'PthH2' 'dHr' 'Pr' 'P0' 'Pth' 'Pthw0' 'eta' 'etaw0'; ...
    '[W]' '[kJ molH2^-1]' '[W]' '[W]' '[W]' '[W]' '[-]' '[-]'};
% Para3 = [PthH2 dHr/1000 Pr P0 Pth Pthw0 Eta Etaw0];
Para3 = [PthH2 dHr/1000 Pr P0 Pth Pthw0 eta etaw0];

Para4Header = {'tdTA' 'tdTE' 'dtdT' 'EdT (>65°C)' 'tSA' 'tSE'; ...
    '[min]' '[min]' '[min]' '[kJ]' '[min]' '[min]'};
Para4 = [tdTAoutput tdTEoutput-mod(tdTEoutput,dt) dtdT EdT ...
    tSAoutput-mod(tSAoutput,dt) tSEoutput-mod(tSEoutput,dt)];

Para5Header = {'Chargennummer' 'Bodenmodul' 'Mittelmodul'}';
Para5 = {CH BM MM}';

xlswrite(outname, Para1Header, 'Ergebnisse', 'A1');
xlswrite(outname, Para1, 'Ergebnisse', 'A3');
xlswrite(outname, Para2Header, 'Ergebnisse', 'A5');
xlswrite(outname, Para2, 'Ergebnisse', 'A7');
xlswrite(outname, Para3Header, 'Ergebnisse', 'A9');
xlswrite(outname, Para3, 'Ergebnisse', 'A11');
xlswrite(outname, Para4Header, 'Ergebnisse', 'A13');
xlswrite(outname, Para4, 'Ergebnisse', 'A15');
xlswrite(outname, Para5Header, 'Ergebnisse', 'J1');
xlswrite(outname, Para5, 'Ergebnisse', 'K1');