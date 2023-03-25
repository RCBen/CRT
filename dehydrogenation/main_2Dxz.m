function [HGout, QH2out, prcalc, Trcalc, dpr] = main_2Dxz(mc, Q, pr, Tr, Tw)
    % [HGout, QH2out, prcalc, Trcalc, dpr] = main_2Dxz(mc, Q, pr, Tr, Tw)
    %
    % hydrogenation grade    @ outlet (-)
    % volume flow rate of H2          (ml min^-1)
    % calculated reaction temperature (K)
    % calculated reaction pressure    (bar)
    % pressure drop in reactor        (Pa)
    % 
    % mc (g)            scalar
    % Q (ml min^-1)     scalar
    % pr (bar)          scalar
    % Tr / Tw (K)       scalar
    %
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 1275 ff)
    
    if nargin == 0
        mc = 200;
        Q = 10;
        pr = 1;
        Tr = 310+273.15;
        Tw = Tr+5;
    end    
    
    mydir = pwd;
    idcs = strfind(mydir,'\');
    p.subdir = mydir(1:idcs(end-3)-1);

    % constants
    if isempty(which('getConstants'))      
        addpath(strcat(p.subdir,'\functions'))
    end
    const = getConstants;  
    
    %% Definition
    p = myDefinitions(p);                                                   % includes p.rate, p.ir, p.irT, p.film, p.iso, p.DPe, p.FPe, p.DEff, p.lambdaPe, p.lambdaEff
    p.plot = 0;                                                             % plot results                  0 off
                                                                            %                               1 on
    p.hw = 3;                                                               % heat transfer coeff. (wall)   1 plate
                                                                            %                               2 packed bed
                                                                            %                               3 bubble column
                                                                            %                               4 Inf (T(z=0) = Tw)
    p.dpr = 1;                                                              % pressure drop coeff. (dpr)    1 Ergun
                                                                            %                               2 Molerus
                                                                            
    if ~isfield(p,'warn') && p.hw == 2 && p.lambdaEff ~= 5
        p.lambdaEff = 5;
        p.warn = fprintf(['p.lambdaEff was set to 5 due to used ' ...
                   'heat transfer correlation (otherwise set p.hw ~= 2)!\n']);
    end
    
    % variables    
    p.db = 0.5e-03;                                                         % bubble diameter (m)
    p.mc = mc;                                                              % mass of catalyst in reactor (g)
    p.dpore = const.dpore;                                                  % average pore diameter (nm)
    
    if mc == 200
        p.HR = 0.0135;                                                          % filling height of reactor (m)
    elseif mc == 400
        p.HR = 0.0270;                                                          % filling height of reactor (m)
    end
    p = myReactor(const, p);
    p = myCoeff(p);
    
    p.nx = 15;                                                            % number of volume elements in x-direction (-)
%     p.nx = ceil(p.LRxz./(2*const.rd));                                      % number of volume elements in x-direction (-)
    p.x = linspace(0, p.LRxz, p.nx);                                        % position of control volumes in x-direction (m)
    p.dx = p.LRxz/p.nx;                                                     % length of volume element in x-direction (equidistant) (m)
    
    p.nz = 2;                                                              % number of volume elements in z-direction (-)  
%     p.nz = ceil(p.HR./(2*const.rd));                                        % number of volume elements in z-direction (-)  
    p.z = linspace(0, p.HR, p.nz);                                          % position of control volumes in z-direction (m)
    p.dz = p.HR/p.nz;                                                       % length of volume element in z-direction (equidistant) (m)
    
    % parameters
    p.HGin = 1;                                                             % hydrogenation grade @ inlet (-)
    p.Q = Q;                                                                % volume flow rate of liquid (ml min^-1)
    p.tau = p.V0/p.Q;                                                       % hydrodynamic residence time (min)
    p.pr = pr;                                                              % reaction pressure (bar)
    p.Tr =  Tr;                                                             % reaction temperature (K)    
    if ~isscalar(Tw)
        idx = ones(size(p.dxF));
        for i = 2:numel(p.dxF)
            idx(i) = find(p.x >= p.dxF(i), 1, 'first');
%             p.Tw(idx(i-1):idx(i)) = Tw(i-1);                                % constant
            p.Tw(idx(i-1):idx(i)) = Tw(i-1) + ...                           % linear interpolation
                (Tw(i)-Tw(i-1))/(idx(i)-idx(i-1)).*(0:(idx(i)-idx(i-1)));
        end
        p.Tw(idx(i):p.nx) = Tw(i);
    else
        p.Tw = repmat(Tw, size(p.x));
    end  
    p.mL0 = p.V0.*(getRho(const, p.pr, p.Tr, p.HGin)/1000);                 % initial mass of liquid phase (g)
    p.nL0 = p.mL0./(p.HGin*const.M18 + (1-p.HGin)*const.M0);                % initial molar amount of liquid phase (mol)

    p.mPt = p.mc.*const.wPtc;                                               % mass of platin (g)
    p.xPt = p.mPt./(p.nL0*const.MPt);                                       % molar fraction of Pt regarding to nL0 (-)  

    %% Calculations
    syms T HG
    rho = getRho(const, [], T, HG);
    p.betaFun = matlabFunction(-1./rho.*diff(rho, T), 'Vars', [T HG]);      % function for thermal expansion coefficient (K^-1)
    
    p.Tamb = const.T0;                                                      % ambient temperature (K)

%     p_estxz = Cooling(p);
    p_estxz = [1.779277560505660e-03, 1.465908578133450e-01];               % fit to experimental data on cooling @ 1 to 3.2 bar (average value)                                                            
    p.kamb = matlabFunction(poly2sym(p_estxz, T));                          % function for ambient heat transfer coefficient (W m^-2 K^-1)
                                                                            % kamb is referred only to reaction volume (without reactor and heating plate)

    % input variables
    m0in = (p.Q/60*1e-06)*getRho(const, [], const.T0, p.HGin);              % mass flow rate in x-direction @ inlet (kg s^-1)
    rhou0in = m0in./p.Ayz;                                                  % mass flux in x-direction @ inlet (kg m^-2 s^-1)       
    w18in = p.HGin*const.M18./(p.HGin*const.M18 + (1-p.HGin)*const.M0);     % mass fraction of H18-DBT @ inlet (-)
    Trin = p.Tr;                                                            % reaction temperature @ inlet (K)
    
    % solve system of ODEs
    x0 = repelem([rhou0in w18in Trin], 1, p.nz);                            % combine initial conditions in one vector
    
    JPattern = [zeros(p.nz, 3*p.nz); ...
        eye(p.nz) ...
        spdiags([ones(p.nz,1) ones(p.nz,1) ones(p.nz,1)], [-1, 0, 1], ...
        p.nz, p.nz) zeros(p.nz); ...
        eye(p.nz) ...
        spdiags([ones(p.nz,1) ones(p.nz,1) ones(p.nz,1)], [-1, 0, 1], ...
        p.nz, p.nz) ...
        spdiags([ones(p.nz,1) ones(p.nz,1) ones(p.nz,1)], [-1, 0, 1], ...
        p.nz, p.nz)];
    if p.iso == 1
        JPattern(2*p.nz+1:end,:) = 0;
    end
    options = odeset('NonNegative', 1:1:3*p.nz, 'JPattern', JPattern);
%     tic;
    [xout, Xout] = ode15s(@odes_2Dxz, p.x, x0, options, const, p);
%     toc;

    % disassemble the solution vector/matrix  
    rhou0 = Xout(:, 0*p.nz+1:1*p.nz)';                                      % mass flux in x-direction (kg m^-2 s^-1)
    w18 = Xout(:, 1*p.nz+1:2*p.nz)';                                        % mass fraction of H18-DBT (-)
    Tr = Xout(:, 2*p.nz+1:3*p.nz)';                                         % reaction temperature (K)
    
    HG = w18./const.M18./(w18./const.M18 + (1-w18)./const.M0);              % hydrogenation grade (-)

    rho = getRho(const, [], Tr, HG);                                        % density of liquid mixture, wo H2 (g l^-1)
    eta = getEta(const, [], Tr, HG)*1e-03;                                  % dynamic viscosity of liquid mixture, wo H2 (Pas)
    u0 = rhou0./rho;                                                        % superficial velocity of liquid in x-direction (m s^-1)
    u0(p.z == 0,:) = 0;                                                     % superficial velocity of liquid in x-direction @ wall (m s^-1)
    
%     Qm = zeros(size(HG));                                                   
%     for x = 1:p.nx
%         Qm(:,x) = myQm(const, p, p.pr, Tr(:,x)', HG(:,x)');                 
%     end
    Qm = myQm(const, p, p.pr, Tr(:)', HG(:)');                              % mass source (g l_bed^-1 s^-1)
    Qm = reshape(Qm, size(Tr));
    mH2 = trapz(xout, trapz(p.z,-Qm))*p.BRxz*1e03;                          % mass flow rate of H2 (g s^-1)
%     mH2 = sum(-Qm,'all')*(p.dx*p.dz)*p.BRxz*1e03;
                   
    nH2max = p.HGin.*(m0in*1e03)*const.wH218/const.MH2;                     % maximum releasable H2 from m0in (mol s^-1)
    HGout = 1-(mH2./const.MH2)/nH2max;                                      % hydrogenation grade @ outlet (-)
    QH2out = (mH2*60*1e03)/const.rhoH2n;                                    % hydrogen volume flow rate @ outlet (ml min^-1)
        
%     rhoH2Id = ((p.pr*1e+05).*const.MH2./(const.R.*Tr))*1e-03;               % density of ideal gas H2 (g l^-1)
%     if p.db <= 2e-03                              
%         phiH2 = 1/2*(1 - sqrt(1 - 4*(12.*eta.* ...                          % volume fraction of H2 (-)
%             cumtrapz(p.z, -Qm))./ ...     
%             (rhoH2Id.*rho.*const.g.*p.db.^2)));
%     end
%     phi = 1-phiH2;                                                          % volume fraction of liquid phase (-)
    
    switch p.dpr
        case 1
            kappa = (150/(2*const.rd)^2*(1-p.epsBed)^2/p.epsBed^3)^-1;
            betaF = 1.75/(2*const.rd).*rho.*(1-p.epsBed)/p.epsBed^3;
            dprdx = -(eta./kappa + betaF.*abs(u0)).*u0;
        case 2
            Rep0 = rho.*(u0./p.epsBed).*(2*const.rd)./eta; 
            r0delta = (0.95/(1-p.epsBed).^(1/3)-1).^-1;
            Eu = 24./Rep0.*(1 + 0.692*(r0delta + 0.5*r0delta.^2)) + ...
                4./sqrt(Rep0).*(1+0.12*r0delta.^1.5) + ...
                (0.4 + 0.891*r0delta.*Rep0.^-0.1);
            dprdx = -3/4*Eu.*rho.*u0.^2./(2*const.rd).* ...
                (p.epsBed.^2./(1-p.epsBed)).^-1;
    end
    dprdx = dprdx - (Qm./p.epsBed.^2).*u0;                                  % pressure drop in x-direction (Pa m^-1)
    dpr = trapz(xout, dprdx')';                                             % pressure drop in reactor (Pa)
    prcalc = (p.pr*1e05 + (cumtrapz(xout, dprdx')' - dpr) + ...             % calculated reaction pressure (bar)
         rho.*const.g.*(p.HR - p.z)')*1e-05;

    Trcalc = mean2(Tr(dsearchn(p.z', p.zE(1)'), dsearchn(xout, p.xE')));    % calculated reaction temperature (K)
    
    if p.plot == 1
        refz = [p.nz, ceil(p.nz/2), 1];
        
        % start 1st figure
        figure()
        grid on
        hold on                                                             % several graphs in one plot
        
        title('spatial behavior of hydrogenation grade');
        xlabel('hydrogenation grade [-]');
        
        plot(HG(:,xout == 0), p.z*1e03,'b-');                               % plot hydrogenation grade @ first x over z-coordinate
        plot(HG(:,xout == p.LRxz), p.z*1e03, 'r-');                         % plot hydrogenation grade of last x over z-coordinate  
        plot(HG(:,2:ceil(p.BRxz/p.dx):end-1), p.z*1e03,'k-');               % plot hydrogenation grade for every x over z-coordinate       
        ylabel('position inside reactor, z [mm]');
        
        legend('x = 0', 'x = L');                                           % draw a legend with the first and last time step
        
        % start 2nd figure
        figure();
        grid on
        
        title(['spatial behavior of superficial velocity ' ...
            'and pressure in reactor']);
        xlim([0 p.LRxz]*1e03)
        xlabel('position inside reactor, x [mm]');

        yyaxis left
        plot(xout*1e03, u0(refz,:)*1e03);
        ylim([0.95*min(u0(:)),1.05*max(u0(:))]*1e03)
        ylabel('superficial velocity of liquid in x-direction [mm/s]');   
        
        yyaxis right
        plot(xout*1e03, (prcalc(refz,:) - p.pr)*1e03);
        ylim([0.95*min(prcalc(:)-p.pr), 1.05*max(prcalc(:)-p.pr)]*1e03)    
        ylabel('gauge pressure [mbarü]'); 

        legend('z = H', 'z = H/2', 'z = 0');                                % draw a legend with the first and last time step
        
        % start 3rd figure
        figure();
        grid on
        hold on
        
        title(['spatial behavior of hydrogenation grade ' ...
            'and reaction temperature in reactor']);
        xlim([0 p.LRxz]*1e03)
        xlabel('position inside reactor, x [mm]');
        
        yyaxis left
        plot(xout*1e03, HG(refz,:));
        ylim([0.95*min(HG(:)), p.HGin])    
        ylabel('hydrogenation grade [-]'); 
        
        yyaxis right
        plot(xout*1e03, Tr(refz,:)-273.15);
        plot(xout*1e03, p.Tw-273.15);
        ylim([0.95*min(Tr(:)-273.15), max(p.Tw-273.15)])    
        ylabel('reaction temperature [°C]');

        legend('z = H', 'z = H/2', 'z = 0','','','', 'Twall');                                % draw a legend with the first and last time step
    end
end