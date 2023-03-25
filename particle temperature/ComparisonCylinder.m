clear variables
% if ~exist('p', 'var') || ~isfield(p,'subdir')
    mydir = pwd;
    idcs = strfind(mydir,'\');
    p.subdir = mydir(1:idcs(end-2)-1);
    addpath(strcat(p.subdir,'\functions'))
% end
const = getConstants;
p = myDefinitions(p);
p = myCoeff(p);

p.ir = 1;
p.irT = 1;
p.film = 1;
p.F = 0;

% p.DPe = 3;                                                                  % correction factor (DPe)       1 ratio 1/10
                                                                            %                               2 Dumanski
                                                                            %                               3 Hugo
                                                                            %                               4 Bruggemann
                                                                            %                               5 Maxwell
                                                                            %                               6 Millington/Quirk     
% p.lambdaPed = 2;                                                             % correction factor (lambdaPe)  1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)
                                                                            %                               3 Harriot
p.FitAS = 0;                                                                % fit to AS data (Diss p.80) in order to evaluate best models
p.deltaMax = 0;                                                             % get deltaMax for each temperature -> delta is divided into more steps


pr = 1;

E = load('EvaluationParticle.mat', 'HG0', 'HG', 'nL0', 'dH', 'dS', 'nd', ...
    'wPtLH', 'wPtLS', 'nTheta', 'mc', 'S', 'H', 'Results', 'Thetarstr');

p.dpore = const.dpore;


if p.FitAS == 1
%     dH = [10 20 140 340 740 1200 3000]*10^-6;
% elseif p. FitAS == 2
    dH = [26 66 95 177]*10^-6;
elseif p.deltaMax == 1
    dH = (193:1:195)*10^-6;
else
    dH = E.dS;
end
phiPt = 1-(const.r./(const.r+dH)).^3;                                       % volume fraction of catalyst with Pt
wPtL = E.wPtLS;                                                             % mass fraction of Pt regarding to mass of active layer for spheres (-)

Vc = E.S.mc/const.rhoBed*1000;                                              % volume of catalyst in bed (thermocouple + spheres) (ml)  

ddH = 5;
dHstart = 15;
dHend = 200;
ndH = (dHend-dHstart)./ddH+1;
dHcalc = linspace(dHstart, dHend, ndH);


% p.epsBed = const.epsBed;
if p.deltaMax == 1
    nL0 = 0.1;
    HG0 = 0.95;
else
    nL0 = E.nL0;
    HG0 = E.HG0;
end
HG = 0.9;                                                                   % allowed values: 0.6:0.05:0.9
X = 1-HG/HG0;
n18 = HG*nL0;
n0 = nL0 - n18;

h = find(abs(E.HG-HG) < 1e-6)+1;
dTexpH = zeros(size(E.Results.S.HG));

Tr = zeros(1,numel(E.Thetarstr));
for j = 1:E.nTheta
    Tr(j) = str2double(E.Thetarstr{j})+273.15;
end

% if p.FitAS == 1
%     Tr = Tr+10;
% end
for i = 1:E.nd
    for j = 1:E.nTheta
        if numel(E.Results.S.dT{i,j}) >= h && ~isempty(E.Results.S.dT{i,j}(h)) 
            dTexpH(i,j) = E.Results.S.dT{i,j}(h);                           % dT for cylinder
        end
    end
end
riS = const.r;
raS = riS + E.dS;
riH = const.r;
raH = riH + E.dH;
wPtS = E.wPtLS;
wPtH = E.wPtLH;
factor = 2*((1./riS-1./raS)./(3*log(raH./riH))).*(((raS.^3-riS.^3))./(raH.^2-riH.^2));
factor2 = wPtS./wPtH;

dTexp = dTexpH.*factor.*factor2;
dTexpU = dTexp + 2*((0.3+0.005*(repmat(Tr, E.nd,1)-273.15))/3);             % dTexp with uncertainties from temperature difference measurement (worst case)

% for l = 1:3
%     p.lambdaPe = l;
    name1 = strcat('lambdaPe', num2str(p.lambdaPe));
%     for k = 1:7
%         p.DPe = k;
        name2 = strcat('DPe', num2str(p.DPe));
%         for i = 1:ndH
        for i = 1:numel(dH)   
            for j = 1:E.nTheta
%                 j = 4;
                p.phiPt = phiPt(i);
                if p.deltaMax == 1
                    % assumption constant cat mass
                    p.Vc = mean(E.S.mc/const.rhoBed*1000);
                    p.wPtL = mean(E.wPtLH);
                    p.wPtc = p.wPtL.*p.phiPt;

                    % assumption constant xPt and wPtc or wPtL
                    p.wPtL = 2.5/100;                                      % goal 2%, but often 2.2 to 2.3%
                    p.wPtc = p.wPtL.*p.phiPt;
%                     p.wPtc = 0.3/100;                                      % between 0.3% and 0.6%
%                     p.wPtL = p.wPtc./p.phiPt;
%                     
                    p.xPt = 0.1/100;
                    p.mc = (p.xPt.*nL0.*195.084)./p.wPtc;
                    p.Vc = p.mc/const.rhoBed*1000;
                else                    
                    p.wPtL = wPtL(i);
                    p.Vc = Vc(i);
                    
%                     p.wPtL = 2.2/100;                                      % goal 2%, but often 2.2 to 2.3%
%                     p.wPtc = p.wPtL.*p.phiPt;
% 
%                     p.xPt = 0.1/100;
%                     p.mc = (p.xPt.*nL0.*195.084)./p.wPtc;
%                     p.Vc = p.mc/const.rhoBed*1000;
                end
                    [reff4, V, etac, P, p] = myReactionRate(const, p, pr, Tr(j), n18, n0, dH(i));
                    R.(name1).(name2).r{i,j} = P(:,1);
                    R.(name1).(name2).R(i,j) = P(end,1);
                    R.(name1).(name2).c18{i,j} = P(:,2);
                    R.(name1).(name2).T{i,j} = P(:,3);
                    R.(name1).(name2).Trs{i,j} = P(end,3);
                    R.(name1).(name2).dTcalc(i,j) = Tr(j)-P(1,3);           % i = diameter, j = temperature, k = DPe model, l = lambdaPe model
                    R.(name1).(name2).etac(i,j) = etac;
                    R.(name1).(name2).reff4(i,j) = reff4;
                    R.(name1).(name2).V(i,j) = V;
            end  
        end
%         R.(name1).(name2).dTcalc(R.(name1).(name2).dTcalc==0) = NaN;
%         if k == 1
%             R.(name1).(name2).dTError = (R.(name1).(name2).dTcalc-dTexp)./dTexp;
%         else
%             R.(name1).(name2).dTError = (R.(name1).(name2).dTcalc-dTexpU)./dTexpU;
%         end
%     end
% end

% save('ComparisonCylinder.mat')

dTcalc = R.(name1).(name2).dTcalc;
dTcalc(dTcalc==0) = NaN;
etac = R.(name1).(name2).etac;
reff4 = R.(name1).(name2).reff4;
V = R.(name1).(name2).V;


%compare Prod
n18 = HG0*nL0;
n0 = nL0 - n18;

wPtc = wPtL.*phiPt;
mPt = E.S.mc.*wPtc;
% Vc = E.S.mc/const.rhoBed*1000;  

p.Tr = Tr;
p.pr = pr;
p.pr = repelem(p.pr,numel(p.Tr));                                   % reaction pressure (bar)  
n18 = repelem(n18,numel(p.Tr));
n0 = repelem(n0,numel(p.Tr));
p.J = numel(n18);

t = linspace(0, 2000, 100);
idx90 = zeros(size(dH,2),size(n18,2));
% idx70 = zeros(size(dH,2),size(n18,2));
dt = zeros(size(dH,2),size(n18,2));
for i = 1:numel(dH)
    p.Vc = Vc(i);
    p.dH = dH(i);
    [~, nrvcalc] = ode15s(@(t, nrv, p) myMolarChange(t, nrv, p), t, ...
        [n18 n0], [], p);
    
    n18rvcalc = nrvcalc(:,1:p.J);
    n0rvcalc = nrvcalc(:,p.J+1:end); 
    X = 1-n18rvcalc./n18rvcalc(1,:);
    HGcalc = n18rvcalc./n18rvcalc(1,:).*HG0;
    for ii = 1:size(HGcalc, 2)
        idx90(i,ii) = find(HGcalc(:,ii) <= 0.9, 1);
%         idx70(i,ii) = find(HGcalc(:,ii) <= 0.7, 1);
%         dt(i,ii) = t(idx70(i,ii)) - t(idx90(i,ii));
        dt(i,ii) = t(idx90(i,ii));
    end
end
if p.deltaMax ~=1
    Prod = (0.95-0.9).*nL0*9.*const.MH2./dt./mPt;
    Prodc = 2.06./Prod(2,1).*Prod;
    disp([round(dH*10^6,-1) Prodc])     
end

if p.FitAS == 1      
        HG0 = 1;
        nL0 = 0.1;
%         p.xPt = 0.1/100;
%         nPt = p.xPt*nL0;
%         mPt = nPt.*const.MPt;
        wPtL = 2.2/100;
        wPtc = wPtL.*phiPt;
%         mc = mPt./wPtc;
%         Vc = mc/const.rhoBed*1000;
        n18 = HG0*nL0;
        n0 = nL0 - n18;

        VPe = 4/3*pi*((const.r+dH).^3)*1000;                                % volume of pellet (l)
        NPe = 50;
        Vc = NPe.*VPe./(1-const.epsBed)*1000;                               % cat volume in ml
        mc = Vc.*const.rhoBed/1000;                                         % cat mass in g
        mPt = mc.*wPtc;
        nPt = mPt./const.MPt;
        p.xPt = nPt./nL0;

        p.Tr = Tr;
        p.pr = pr;
        p.pr = repelem(p.pr,numel(p.Tr));                                   % reaction pressure (bar)  
        n18 = repelem(n18,numel(p.Tr));
        n0 = repelem(n0,numel(p.Tr));
        p.J = numel(n18);
        
        t = linspace(0, 6000, 200);
        idx90 = zeros(size(dH,2),size(n18,2));
        idx70 = zeros(size(dH,2),size(n18,2));
        dt = zeros(size(dH,2),size(n18,2));
        for i = 1:numel(dH)
            p.Vc = Vc(i);
            p.dH = dH(i);
            [~, nrvcalc] = ode15s(@(t, nrv, p) myMolarChange(t, nrv, p), t, ...
                [n18 n0], [], p);
            
            n18rvcalc = nrvcalc(:,1:p.J);
            n0rvcalc = nrvcalc(:,p.J+1:end); 
            X = 1-n18rvcalc./n18rvcalc(1,:);
            HGcalc = n18rvcalc./n18rvcalc(1,:).*HG0;
            for ii = 1:size(HGcalc, 2)
                idx90(i,ii) = find(HGcalc(:,ii) <= 0.9, 1);
                idx70(i,ii) = find(HGcalc(:,ii) <= 0.7, 1);
                dt(i,ii) = t(idx70(i,ii)) - t(idx90(i,ii));
            end
        end
        
        Prod = (0.9-0.7).*nL0*9.*const.MH2./dt./mPt;
        Prodc = 2.4./Prod(1,1).*Prod;
        disp(Prodc)        
 
%     Prod = etac.*reff4.*V;
%     Prodc = 7.2./Prod(1,1).*Prod(:,1);
%     disp(flip(Prodc',2))
% elseif p.FitAS == 2
%     mPtS = 0.0789.*dH'*1e6-0.0748;                                           % fit to data mPt vs d for xPt in layer approx. 2 to 2.2
%     mPtS = 0.0713.*dH'*1e6+1.3195;
%     wPtc = 2.2/100.*phiPt';
%     Prod = etac.*reff4.*V./repmat(mPtS, 1, numel(Tr));
%     Prod = etac.*reff4.*V./repmat(wPtc, 1, numel(Tr));                            % mc is constant, thus wPtc can be used instead
%     Prodc = 2.4./Prod(1,1).*Prod;
%     disp(Prodc)
elseif p.deltaMax == 1
    disp([dH'*10^6 etac<0.95]);
end

%% plot
xmin = 10;
xmax = 210;
ymin = 0;
if p.FitAS ~=1
    ymax = 2;
else
    ymax = 15;
end
ylabelname = 'dT';
ylabelunit = '(K)';
legendnameD = {'20 um', '40 um', '60 um', '80 um', '110 um', '150 um', '200 um'};
legendnameT = {'310 °C','300 °C', '290 °C', '280 °C'};

% if p.DPe == 1
   dTexpPlot = dTexp;
% else
%    dTexpPlot = dTexpU; 
%     ymax = 3.5;
% end

if p.deltaMax == 0
    figure()
    
    title ('Calculated and experimental results');
    %C={'[.88 .41 .12]';'[.85 .35 .11]';'[.82 .27 .11]';'[.79 .19 .11]'};
    semilogx((R.(name1).(name2).R-const.r)*10^6, dTcalc(:,1),'--','LineWidth',1.5)
    hold on
    semilogx((R.(name1).(name2).R-const.r)*10^6, dTcalc(:,2),'--','LineWidth',1.5)
    semilogx((R.(name1).(name2).R-const.r)*10^6, dTcalc(:,3),'--','LineWidth',1.5)
    semilogx((R.(name1).(name2).R-const.r)*10^6, dTcalc(:,4),'--','color','[0 0.55 0]','LineWidth',1.5)
    set(gca,'ColorOrderIndex',1);
    semilogx(E.dS*10^6, dTexpPlot(:,1),'x','LineWidth',1.5)
    semilogx(E.dS*10^6, dTexpPlot(:,2),'d','LineWidth',1.5)
    semilogx(E.dS*10^6, dTexpPlot(:,3),'o','LineWidth',1.5)
    semilogx(E.dS*10^6, dTexpPlot(:,4),'v','color','[0 0.55 0]','LineWidth',1.5)
    xlim([xmin xmax])
    ylim([ymin ymax])
    xticks([20 40 60 80 100 200])
    set(gca,'FontSize',14)
    
    xx =xlabel('Porous layer thickness [µm]');
    yy =ylabel('Temperature drop [C°]');
    h = legend('T = 310°C','T = 300°C','T = 290°C','T = 280°C','Location','northwest');
    set(h,'Fontsize',12)
    set(xx,'Fontsize',14)
    set(yy,'Fontsize',14)
    
    grid on
    hold off
end

%---------------------------------------------------
%How to plot Temperature drop vs. changing d at T=310 (when i=4 in {j,i})
if p.deltaMax == 0 && p.FitAS == 0
    figure();
    % C={'[.97 .57 .12]';'[.92 .47 .12]';'[.88 .41 .12]';'[.85 .35 .11]';'[.82 .27 .11]';'[.79 .19 .11]'};
    j = 1;
    for i = 1:E.nd
        plot ((R.(name1).(name2).r{i,j}-const.r)*10^6, R.(name1).(name2).T{i,j}-273.15, '-d','LineWidth',1.5);
        hold on
    end
    xlim([xmin xmax])
    ylabel('Temperature') % left y-axis 
    xlabel('Porous layer (mm)');
    axis tight;
    grid on
    %---------------------------------------------------
    %How to plot Concentration drop vs. changing d at T=310 (when i=4 in {j,i})
    figure();
    C={'[.97 .57 .12]';'[.92 .47 .12]';'[.88 .41 .12]';'[.85 .35 .11]';'[.82 .27 .11]';'[.79 .19 .11]'};
    for i = 1:E.nd
        plot ((R.(name1).(name2).r{i,j}-const.r)*10^6, R.(name1).(name2).c18{i,j}, '-x','LineWidth',1.5);
        hold on
    end
    xlim([xmin xmax])
    ylabel('Concentration') 
    xlabel('Porous layer (um)');
    axis tight;
    grid on
end




%%% nested function

function dndt = myMolarChange(~, nrv, p)       
%     p.J = numel(nrv)/2;
    n18 = nrv(1:p.J)';                                                      % molar amount of H18-DBT in liquid phase (mol)
    n0 = nrv(p.J+1:end)';                                                   % molar amount of H18-DBT in liquid phase (mol)  
    
    const = getConstants;
    [reff4, V, ~, ~, p] = myReactionRate(const, p, p.pr, ...                % effective reaction rate (mol l_bed^-1 s^-1) & volume of liquid phase (l)
        p.Tr, n18, n0, p.dH);       
    dn18dt = (const.nu(1).*(V+p.Vc.*(1-const.epsBed)/1000).*reff4)';        % molar change of H18-DBT in liquid phase (mol s^-1)
    dn0dt = (const.nu(2).*(V+p.Vc.*(1-const.epsBed)/1000).*reff4)';         % molar change of H0-DBT in liquid phase (mol s^-1)

%     [reff4, ~, ~, ~, p] = myReactionRate(const, p, p.pr, ...                % effective reaction rate (mol l_bed^-1 s^-1) & volume of liquid phase (l)
%         p.Tr, n18, n0, p.dH);       
%     dn18dt = (const.nu(1).*(p.Vc/1000).*reff4)';                            % molar change of H18-DBT in liquid phase (mol s^-1)
%     dn0dt = (const.nu(2).*(p.Vc/1000).*reff4)';                             % molar change of H0-DBT in liquid phase (mol s^-1)

    dndt = [dn18dt; dn0dt];
end


