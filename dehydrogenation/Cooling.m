function [p_estxz, Tcalc, kambxz] = Cooling(p)
    % Bi_max = 0.0042 -> lambdaEff, i.e. dT in fluid, can be neglected
    %
    
    if isempty(which('getConstants'))
        mydir = pwd;
        idcs = strfind(mydir,'\');
        p.subdir = mydir(1:idcs(end-2)-1);
        addpath(strcat(p.subdir,'\functions'))
    end
    
    if nargin == 0
        p.DEff = 4;                                                         % correction factor (DEff)      1 ratio 1/10
                                                                            %                               2 Dumanski
                                                                            %                               3 Hugo
                                                                            %                               4 Bruggemann
                                                                            %                               5 Maxwell
                                                                            %                               6 Millington/Quirk 
        p.rhocpPed = 2;                                                     % correction factor (rhocpPed)  1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)
                                                                            %                               3 Harriot                                                                    
        p.rhocpPe = 1;                                                      % correction factor (rhocpPe)   1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)
        p.rhocpEff = 1;                                                     % correction factor (rhocpEff)  1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)
                                                                            %                               3 power law
                                                                            %                               4 Zehner/Bauer/Schlünder (primary)
                                                                            %                               5 Zehner/Bauer/Schlünder (secondary)
    else
        p.rhocpPed = p.lambdaPed;
        p.rhocpPe = p.lambdaPe;
        p.rhocpEff = p.lambdaEff;
    end
    
    %% constants
    const = getConstants;
    
    %% variables    
    if nargin == 0
        p.mc = 200; 
        p.Tamb = const.T0;
        p = myReactor(const, p);
    end
    
    samples = [1, 2, 3];
    nsamples = numel(samples);

    %% experimental data (E6 - E10)
    HG = [1, 0.8102, 0.7540];
    HG = HG(samples);
    t0 = 0;
    tend = 55800;
    dt = 1800;
    t = linspace(t0, tend, 55800/dt + 1);
    
    pr = [1, 2.7, 3.2];
    pr = pr(samples);
    
    % pr = 1 bar
    A005 = [267.7 250.98 235.38 221.52 209.44 198.14 187.8 177.86 168.7 ...
        160.06 151.68 144.36 137.06 130.48 124.06 118.16 112.84 107.66 ...
        102.82 98.4 94.22 90.24 86.46 83.02 79.74 76.7 73.88 71.16 68.6 ...
        66.18 63.92 61.76];
    
    % pr = 2.7 bar
    A009 = [272.56 256.18 241.2 227.44 214.66 202.82 191.78 181.52 171.94 ...   
        162.66 154.14 146.2 138.8 131.82 125.34 119.28 113.62 108.32 ...
        103.38 98.74 94.42 90.38 86.58 83.06 79.7 76.6 73.64 70.84 68.3 ...
        65.82 63.52 61.38];
    
    % pr = 3.2 bar
    A010 = [272.64 256.24 241.38 227.66 214.94 203.1 192.08 181.8 171.92 ...    
        162.9 154.44 146.48 139.04 132.06 125.58 119.48 113.84 108.54 ...
        103.56 98.94 94.62 90.6 86.78 83.22 79.88 76.74 73.84 71.04 68.42 ...
        66.02 63.72 61.52];

    Tr = [A005; A009; A010] + 273.15;
    Tr = Tr(samples,:);

    p_lb = [];
    p_start = [2 -1];
    p_ub = [];
    D = diag([1e-03, 1]);

    optionsLSQ = optimoptions('lsqnonlin', 'Algorithm', ...
        'trust-region-reflective', 'Display', 'off', ...
        'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', ...
        200*numel(p_start));
    [p_estxz, ~, Residuals] = ...
        lsqnonlin(@(p_est) myModel(p_est*D), p_start, ...
        p_lb, p_ub, optionsLSQ);

    p_estxz = p_estxz*D;
    Tcalc = Tr - Residuals';
    kambxz = polyval(p_estxz, Tcalc);   
    
    %% nested functions
    function Residuals = myModel(p_est)
        options = odeset('NonNegative', size(Tr,1));

        [~, Tcalc] = ode15s(@(t, T, p_est) mydT(t, T, p_est), ...
            t', Tr(:,t==0), options, p_est);
        
        Residuals = Tr' - Tcalc;

        function dTdt = mydT(~, T, p_est)
            T = T';
            kambxz = polyval(p_est, T);
            
            %% compare to getLambdaPe.m
            rho = getRho(const, pr, T, HG);                                 % density of liquid mixture(kg m^-3)
            cp = getCp(const, pr, T, HG)*1e06./ ...                         % isobaric heat capacity of liquid mixture (J kg^-1 K^-1)
                (HG.*const.M18+(1-HG).*const.M0)';       
            [rhoAlox, cpAlox] = getAlox(const, T);                          % density % isobaric heat capacity of alumina (kg m^-3 / J kg^-1 K^-1)
            rhoPt = 21557.19 - 0.5675783*T - 1.7525e-05*T.^2 - ....
                3.171806e-08*T.^3 + 4.698968e-12*T.^4;
            cpPt = 122.2187 + 0.03986346*T - 1.836174e-05*T.^2 + ...
                7.556773e-09*T.^3;
    
            rhocp = rho.*cp;
            rhocp = rhocp((1:nsamples) + ((1:nsamples)-1)*(nsamples+size(rho,2)));

            rhocpAlox = rhoAlox.*cpAlox;
            rhocpPt = rhoPt.*cpPt;
          
            if isfield(p,'phiPt')                                           % volume fraction of catalyst with Pt
                phiPt = p.phiPt;    
            else
                phiPt = const.phiPt;
            end
            if isfield(p,'wPtL')                                            % mass fraction of Pt regarding to mass of active layer     
                wPtL = p.wPtL;          
            elseif isfield(p,'wPtc')
                wPtL = p.wPtc./phiPt;
            else
                wPtL = const.wPtc./phiPt;
            end
            
            phiPtL = wPtL./rhoPt./(wPtL./rhoPt + (1-wPtL)./rhoAlox);        % volume fraction of Pt regarding to volume of active layer  
    
            rhocpL = (1-phiPtL).*rhocpAlox + phiPtL.*rhocpPt;               % rhocp of active layer (J m^-3 K^-1) 
            rhocpL = repmat(rhocpL,size(HG,1),1);
    
            switch p.rhocpPed
                case 1
                    %1 - Parallel slabs / rule of mixtures (upper bound)
                    rhocpPed = rhocpL.*(1-const.epsPe) + const.epsPe.*rhocp;
                case 2
                    %2 - Serial slabs (lower bound)
                    rhocpPed = (const.epsPe./rhocp + (1-const.epsPe)./rhocpL).^-1;
                case 3
                    %3 - Harriot
                    delta = 2.*(0.64*3.*(HG.*const.M18 + ...                % mean free path of liquid = mean molecular spacing (m)
                        (1-HG).*const.M0)/1000./...
                        (4*pi.*getRho(const, pr, T, HG).*const.NA)).^(1/3)*1e09; 
                    d = 0.1*(p.dpore*1e+09);
                    Kn = repmat(delta(samples + (samples-1)*nsamples),1,2)./d;
                    Delta = d./(1+d);
                    chi = (rhocp./(1+2*Kn))./rhocpL;
                    rhocpPed = rhocpL.*(Delta.^2 + chi.*(1-Delta).^2 + ...
                        (1 + d + 1./(chi.*Delta)).^-1);
            end
    
            rhocpPer = repmat(rhocpAlox,size(HG,1),1);                      % effective rhocp of support material (J m^-3 K^-1)    
            switch p.rhocpPe                                                % effective rhocp of porous pellet (catalyst) (J m^-3 K^-1) 
                case 1
                    %1 - Parallel slabs / rule of mixtures (upper bound)
                    rhocpPe = rhocpPer.*(1-phiPt) + phiPt.*rhocpPed;  
                case 2
                    %2 - Serial slabs (lower bound)      
                    rhocpPe = ((1-phiPt)./rhocpPer + phiPt./rhocpPed).^-1;  
            end
    
            rhocpEff = p.epsBed.*rhocp + (1-p.epsBed).*rhocpPe;
            %%
            
            dTdt = (-kambxz*p.SvRxz./rhocpEff.*(T-p.Tamb))';
        end
    end
end