function [mb, p] = myMassBalance(p)
    % [mb, p] = myMassBalance(p)
    % 
    % solve mass balance
    
    %% Definitions
    optionsFsolve = optimset('Display','off');
    
    % constants
    if ~isfield(p,'subdir')
        mydir = pwd;
        idcs = strfind(mydir,'\');
        p.subdir = mydir(1:idcs(end-2)-1);
    end
    if isempty(which('getConstants'))
        addpath(strcat(p.subdir,'\functions'))
    end
    const = getConstants;
    
    % parameters
    p.nP = numel(p.pr);                                                     % number of different absolute pressures (-)
    p.nTheta = numel(p.Tr);                                                 % number of different temperatures (-)
    p.I = numel(p.t);                                                       % number of considered reaction times (-)
    p.J = p.nP*p.nTheta;                                                    % number of considered reactions (-)
    
    % reactor                                                              
    p.pr = repelem(p.pr,p.nTheta);                                          % reaction pressure (bar)                                                            
    p.Tr = repmat(p.Tr,1,p.nP);                                             % reaction temperature (K)
    p.Vr = repelem(p.Vr,p.nTheta);                                          % reactor volume (ml)
    
    % LOHC                                 
    p.nL0 = repelem(p.nL0,p.nTheta);                                        % initial molar amount of liquid phase (mol)
    p.m0 = p.nL0.*(p.HG(1,:)*const.M18 + (1-p.HG(1,:))*const.M0);           % initial mass of liquid phase (g)
    p.V0 = p.m0./(getRho(const, p.pr, p.Tr, p.HG(1,:))/1000);               % initial volume of liquid phase (ml)
                   
    
    % catalyst
    if ~isfield(p, 'mc')
        p.xPt = repelem(p.xPt,p.nTheta);                                    % molar fraction of Pt regarding to nL0 (-)
        p.mPt = p.xPt.*p.nL0*const.MPt;                                     % mass of platin (g)
        if ~isfield(p,'wPtc')                                               % mass fraction of Pt regarding to mass of cat (-) 
            if isfield(p,'wPtL') 
                if isfield(p,'phiPt')                                       % volume fraction of catalyst with Pt
                    phiPt = p.phiPt;    
                else
                    phiPt = const.phiPt;
                end
                wPtc = p.wPtL.*phiPt;
            else
                wPtc = const.wPtc;
            end  
        else 
            wPtc = p.wPtc;
        end
        p.mc = p.mPt/wPtc;                                                   % mass of cat (g)
    end
    
    %% Data import
    p.dt = [0; p.t(2:end) - p.t(1:end-1,:)];                                % time steps (s)    

    %% Calculation
    % assumptions: ideal gas (Dalton)    
    %              perfectly mixed gas (vapour) phase
    %              VLE
    %              inert gases (I = Ar, He, N2) are neglected 
    %                  -> xI = 0
    %                  -> pr = p18 + p0 + pH2 (@p.tr > 0),
    %              xH2 -> 0 -> Henrys law
    %                  -> H2 is supercritical fluid @Tr
    %                  -> pLVH2*phiH2LV = fPH2 = Prausnitz fugacity
    %              perfectly mixed liquid phase & time for sample negligible
    %                  -> x18s(ample) = x18, x0s = x0, xH2s = xH2
    %              sample contains only liquid phase (n18vs = n0vs = 0)
    %              incompressible (small liquid bulk modulus)
    %                  -> rho does not depend on pressure      
    %              HG = x18/(1-xH2) = n18/(n18+n0) (@T0 -> nH2=0 -> HG=x18)   
    %              NO loss of DBT due to evaporation (reflux condenser)
    %              dn18v << dn18r & dn0v << dn0r
    %                  -> x18 = x18r, x0 = x0r
    %              NO entrainment
    %              NO heating (HG(p.tr=0) = 1), but pLV(p.tr=0) > 0: 
    %                  - @=1 bar: cat added @p.tr = 0 (=Tr reached)
    %                  - @>1 bar: heating with cat @pr = 70 bar
    %                             p.tr=0 (=pr reached)    
    
    ZG = ones(p.I, p.J);                                                    % compressibility factor of real gas mixture (-)
    y18 = zeros(p.I, p.J);                                                  % molar fraction of H18-DBT in gas (vapour) phase (-)     
    y0 = zeros(p.I, p.J);                                                   % molar fraction of H0-DBT in gas (vapour) phase (-)
    x18 = p.HG;                                                             % molar fraction of H18-DBT in liquid phase (-)  
    x0 = 1-x18;                                                             % molar fraction of H0-DBT in liquid phase (-) 
    xH2 = zeros(p.I, p.J);                                                  % molar fraction of H2 in liquid phase (-)
    
%     if p.evap == 1
%         vle = VLE(const, p.pr, p.Tr, p.HG, p.real);
%         ZG = vle.ZG;                                                            
%         y18 = vle.y18;                                                          
%         y0 = vle.y0;                                                            
%         x18 = vle.x18;                                                          
%         x0 = vle.x0;                                                            
%         xH2 = vle.xH2;
%     end
    rho = getRho(const, p.pr, p.Tr, x18, x0);                               % density of liquid mixture, containing H2 (g l^-1)    

    % mass, volume and molar amount
    p.Vc = p.mc/const.rhoBed*1000;                                          % volume of catalyst in bed (ml)                                                  
%     p.epsBed = const.epsBed;
    p.NPe = round((1-const.epsBed).*(p.Vc/1000)./const.V);                  % number of pellets in bed (-)
    
    n18 = zeros(p.I,p.J);                                                   % molar amount of H18-DBT in liquid phase (mol)
    n0 = zeros(p.I,p.J);                                                    % molar amount of H0-DBT in liquid phase (mol)
    nH2 = zeros(p.I,p.J);                                                   % molar amount of H2 in liquid phase (mol)
    
    V = zeros(p.I,p.J);                                                     % volume of liquid phase (ml)
    VG = zeros(p.I,p.J);                                                    % volume of gas (vapour) phase (ml)
  
    % evaporation
    dn18v = zeros(p.I,p.J);                                                 % molar change of H18-DBT in gas (vapour) phase due to evaporation (mol)
    dn0v = zeros(p.I,p.J);                                                  % molar change of H0-DBT in gas (vapour) phase due to evaporation (mol)
   
    n18v = zeros(p.I,p.J);                                                  % molar amount of H18-DBT in gas (vapour) phase (mol)
    n0v = zeros(p.I,p.J);                                                   % molar amount of H0-DBT in gas (vapour) phase (mol)    

    % reaction
    dn18r = zeros(p.I,p.J);                                                 % molar change of H18-DBT in liquid phase due to reaction (mol)         
    dn0r = zeros(p.I,p.J);                                                  % molar change of H0-DBT in liquid phase due to reaction (mol)
    dnH2gesr = zeros(p.I,p.J);                                              % molar change of H2 due to reaction (mol)
    
    n18r = zeros(p.I,p.J);                                                  % molar amount of H18-DBT in liquid phase neglecting sample + evaporation (mol)
    n0r = zeros(p.I,p.J);                                                   % molar amount of H0-DBT in liquid phase neglecting sample + evaporation (mol)
    nH2gesr = zeros(p.I,p.J);                                               % molar amount of H2 neglecting evaporation (mol)   

    % fsolve
    exitflag = zeros(p.I,p.J);
    
    % VG = Vr - 1000*((n18*M18 + n0*M0 + nH2*MH2)/rho) - Vc*(1-epsBed);
    % n18v = y18*(pr*1e+05)*(VG*1e-06)/(ZG*R*T);
    % n0v = y0*(pr*1e+05)*(VG*1e-06)/(ZG*R*T);
    % n0 = n00 - n0v;
    % n18 = n180 - n18v;
    % @p.tr = 0: x0 = xH2 = 0 -> nH2 = 0
    HG0 = p.HG(1,:);
    n18r(1,:) = HG0.*p.nL0;
    n0r(1,:) = (1-HG0).*p.nL0;

    VG(1,:) = (p.Vc.*(1-const.epsBed) - p.Vr + ...
        (1000*(n18r(1,:)*const.M18 + n0r(1,:)*const.M0))./rho(1,:))./ ...
        ((1000*(const.M18*y18(1,:) + const.M0*y0(1,:)).*...
        (p.pr*1e+05)*1e-06)./(ZG(1,:)*const.R.*p.Tr.*rho(1,:)) - 1);
    
    n18v(1,:) = y18(1,:).*(p.pr*1e+05).*(VG(1,:)*1e-06)./ ...
        (ZG(1,:)*const.R.*p.Tr);
    dn18v(1,:) = n18v(1,:);
    n18(1,:) = n18r(1,:) - dn18v(1,:);
    
    n0v(1,:) = y0(1,:).*(p.pr*1e+05).*(VG(1,:)*1e-06)./ ...
        (ZG(1,:)*const.R.*p.Tr);
    dn0v(1,:) = n0v(1,:);
    n0(1,:) = n0r(1,:) - dn0v(1,:);
    
    V(1,:) = 1000*((n18(1,:)*const.M18 + n0(1,:)*const.M0)./rho(1,:));

    for i = 2:p.I
        %% reaction + evaporation
        Y0 = [n18(i-1,:); n0(i-1,:); nH2(i-1,:); V(i-1,:); VG(i-1,:); ...
            dn18v(i-1,:); dn0v(i-1,:); n18v(i-1,:); n0v(i-1,:); ...
            dn18r(i-1,:) - dn18v(i-1,:); dn0r(i-1,:) - dn0v(i-1,:)];
        
        for j = 1:p.J
            if j == 1
                Y = zeros(size(Y0,1),p.J);
            end
            [Y(:,j), ~, exitflag(i,j)] = ...
                fsolve(@(x) myMassBalanceFun(x, const, p), Y0(:,j), ...
                optionsFsolve);
        end
        
        n18(i,:) = Y(1,:);
        n0(i,:) = Y(2,:);        
        nH2(i,:) = Y(3,:);
        
        V(i,:) = Y(4,:);
        VG(i,:) = Y(5,:);
        
        dn18v(i,:) = Y(6,:);
        dn0v(i,:) = Y(7,:);
        
        n18v(i,:) = Y(8,:);
        n0v(i,:) = Y(9,:);     
        
        dn18r(i,:) = Y(10,:) + dn18v(i,:);    
        dn0r(i,:) = -dn18r(i,:);
        dnH2gesr(i,:) = 9*dn0r(i,:);       

        n18r(i,:) = n18(i,:) + dn18v(i,:);
        n0r(i,:) = n0(i,:) + dn0v(i,:);
        nH2gesr(i,:) = nH2gesr(i-1,:) + dnH2gesr(i,:);
    end

    p.epsBed = V./(V+p.Vc);
    mges = (n18 + n18v)*const.M18 + (n0 + n0v)*const.M0 + ...               % total mass (for checking conservation of mass) (g)        
        nH2gesr*const.MH2;
    rm = const.nu(1)*dn18r./(p.mc.*p.dt);                                   % reaction rate (mol g^-1 s^-1)
    rm(isnan(rm)) = 0;
    
    %% Saving variables to structure array
    mb.n18 = n18;                                                           % molar amount of H18-DBT in liquid phase (mol)
    mb.n0 = n0;                                                             % molar amount of H0-DBT in liquid phase (mol)
    mb.nH2 = nH2;                                                           % molar amount of H2 in liquid phase (mol)

    mb.x18 = x18;                                                           % molar fraction of H18-DBT in liquid phase (-)                                                          
    mb.x0 = x0;                                                             % molar fraction of H0-DBT in liquid phase (-)    
    mb.xH2 = xH2;                                                           % molar fraction of H2 in liquid phase (-)
    
    mb.rho = rho;
    mb.V = V;                                                               % volume of liquid phase (ml)
    mb.VG = VG;                                                             % volume of gas (vapour) phase (ml)
      
    % evaporation
    mb.dn18v = dn18v;                                                       % molar change of H18-DBT in gas (vapour) phase due to evaporation (mol)
    mb.dn0v = dn0v;                                                         % molar change of H0-DBT in gas (vapour) phase due to evaporation (mol)
   
    mb.n18v = n18v;                                                         % molar amount of H18-DBT in gas (vapour) phase (mol)
    mb.n0v = n0v;                                                           % molar amount of H0-DBT in gas (vapour) phase (mol)    
    
    % reaction
    mb.dn18r = dn18r;                                                       % molar change of H18-DBT in liquid phase due to reaction (mol) 
    mb.dn0r = dn0r;                                                         % molar change of H0-DBT in liquid phase due to reaction (mol)
    mb.dnH2gesr = dnH2gesr;                                                 % molar change of H2 due to reaction (mol)
    
    mb.n18r = n18r;                                                         % molar amount of H18-DBT in liquid phase neglecting sample + evaporation (mol)
    mb.n0r = n0r;                                                           % molar amount of H0-DBT in liquid phase neglecting sample + evaporation (mol)
    mb.nH2gesr = nH2gesr;                                                   % molar amount of H2 neglecting evaporation (mol)   
    
    mb.n18rv = mb.n18r - mb.dn18v;                                          % molar amount of H18-DBT in liquid phase neglecting sample (mol)
    mb.n0rv = mb.n0r - mb.dn0v;                                             % molar amount of H0-DBT in liquid phase neglecting sample (mol)
    
    mb.mges = mges;                                                         % total mass (for checking conservation of mass) (g)
    mb.rm = rm;                                                             % reaction rate (mol g^-1 s^-1)
    
    % fsolve
    mb.exitflag = exitflag;
 
    %% Plot
    plot1Var = sprintf('rm');
    plot1VarUnit = sprintf('mol g^{-1} s^{-1}');
    if nargin == 0
        ymin = sign(min(mb.(plot1Var)(:)))*abs(min(mb.(plot1Var)(:)));
        ymax = sign(max(mb.(plot1Var)(:)))*abs(max(mb.(plot1Var)(:)));
        legendnameP = {'5 bar', '4 bar', '3 bar', '2 bar', '1 bar'};
        legendnameT = {'270 °C','280 °C', '290 °C', '300 °C', '310 °C'};
        xlabelname = 't (s)';
        ylabelname1 = sprintf('%s (%s)', plot1Var, plot1VarUnit);
        
        figure()
        for j = 1:p.nP        
            if p.nP > 1
                subplot(2,round(p.nP/2),j)
            end
            plot(p.t,mb.(plot1Var)(:,linspace(j*p.nTheta, ...
                p.nTheta*(j-1)+1, p.nTheta)),'-');
            title(sprintf('\\bf %d bar',j))
            ylim([ymin,ymax])

            if j == 1
                legend(legendnameT((6-p.nTheta):5))
                ylabel(ylabelname1)
                if p.nP == 1
                    xlabel(xlabelname)                
                end
            elseif j == 2 && p.nP <= 3
                xlabel(xlabelname)           
                if p.nP <= 2
                    ylabel(ylabelname1)                
                end
            elseif j == 3
                xlabel(xlabelname)
                if p.nP <= 4
                     ylabel(ylabelname1)                
                end
            elseif j == 4
                xlabel(xlabelname)
                if p.nP >= 5
                    ylabel(ylabelname1)
                end
            elseif j == 5
                xlabel(xlabelname)
            end
        end
    
        figure()
        for j = 1:p.nTheta        
            if p.nTheta > 1
                subplot(2,round(p.nTheta/2),j)
            end
            plot(p.t,mb.(plot1Var)(:,linspace(p.nTheta*(p.nP-1)+j, j, ...
                p.nP)),'-');
            title(sprintf('\\bf %d °C',p.Tr(j) - 273.15))
            ylim([ymin,ymax])

            if j == 1
                legend(legendnameP((6-p.nP):5))
                ylabel(ylabelname1)
                if p.nTheta == 1
                    xlabel(xlabelname)                
                end   
            elseif j == 2 && p.nTheta <= 3
                xlabel(xlabelname)           
                if p.nTheta <= 2
                    ylabel(ylabelname1)                
                end
            elseif j == 3
                xlabel(xlabelname)
                if p.nTheta <= 4
                     ylabel(ylabelname1)                
                end
            elseif j == 4
                xlabel(xlabelname)
                if p.nTheta >= 5
                    ylabel(ylabelname1)
                end
            elseif j == 5
                xlabel(xlabelname)
            end
        end
    end
    
    %% nested functions
    function Y = myMassBalanceFun(X, const, p)
        Y = [X(1) - (x18(i,j)*(X(1) + X(2) + X(3))); ...                    % n18 = x18*(n18 + n0 + nH2);
            X(2) - n0(i-1,j) - X(11); ...                                   % n0 = n00 + dn0;
            X(3) - (xH2(i,j)*(X(1) + X(2) + X(3))); ...                     % nH2 = xH2*(n18 + n0 + nH2);
            X(4)*rho(i,j) - ...                                             % V = 1000*(n18*M18 + n0*M0 + nH2*MH2)/rho;
                1000*(X(1)*const.M18 + X(2)*const.M0 + X(3)*const.MH2); ...
            X(5) - p.Vr(1,j) + X(4) + p.Vc(1,j).*(1-const.epsBed); ...      % VG = Vr - V - Vc*(1-epsBed);    
            X(6) + n18v(i-1,j) - X(8); ...                                  % dn18v = n18v - n18v0;
        	X(7) + n0v(i-1,j) - X(9); ...                                     % dn0v = n0v - n0v0;    
            X(8)*(ZG(i,j)*const.R*p.Tr(1,j)) - ...                          % n18v = y18*(pr*1e+05)*(VG*1e-06)/(Z*R*Tr);
                (y18(i,j)*(p.pr(1,j)*1e+05)*(X(5)*1e-06)); ...         
            X(9)*(ZG(i,j)*const.R*p.Tr(1,j)) - ...                          % n0v = y0*(pr*1e+05)*(VG*1e-06)/(Z*R*Tr);    
                (y0(i,j)*(p.pr(1,j)*1e+05)*(X(5)*1e-06)); ...
            X(1) - n18(i-1,j) - X(10); ...                                  % n18 = n180 + dn18;     
        	X(11) + (X(10) + X(6)) + X(7)];                                 % dn0 = -(dn18 + dn18v) - dn0v;	
    end
end