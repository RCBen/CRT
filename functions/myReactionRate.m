function [reff4, V, etac, P, p] = myReactionRate(const, p, pr, Tr, n18, n0, d)
    % [reff4, V, etac, P, p] = myReactionRate(const, p, pr, Tr, n18, n0, d)
    %
    % effective reaction rate for dehydrogenation of H18-DBT (mol l_bed^-1 s^-1)
    % volume of liquid phase (l)
    % effectivenes factor of catalyst (-)
    % concentration and temperature profile inside catalyt (mol l^-1) and (K)
    % a priori criteria for heat and mass transport limitation in film (p.MearsH, p.MearsM) and pellet (p.AndersonH) (-)
    % 
    % p.rate (dependency)       0 partial density (rho)
    %                           1 molar concentration (c)
    % p.i(nternal) r(esistance) 0 off
    %                           1 on
    % p.irT (ir regarding T)    0 off
    %                           1 on                                                                        
    % p.film (film limitation)  0 off
    %                           1 on
    % 
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % n18 / n0 (mol)    scalar / vector / matrix
    % d (m)             scalar
 
    p.BiH = 25;                                                             % (heat transfer) Biot number BiH = alphaPe*dPe/lambdaPe (-)
    p.BiM = 1e+12;                                                          % mass transfer Biot number BiM = betaPe*dPe/DPe (-)
                                                                            % set Bim to any arbitrary, very high value (e.g. 1e12) to neglect mass transport limitation in film (if p.film = 1, see myDefinitions.m)
    nP = numel(pr);
    nTheta = numel(Tr);
    if ~isempty(pr)        
        if ~isrow(pr)
            error('pr has to be a scalar or a row-vector')
        elseif isrow(pr) && (nP ~= nTheta || ~isscalar(pr) && pr(1) ~= pr(2))
            pr = repelem(pr,1,nTheta);
            Tr = repmat(Tr,1,nP);
        end
    end
    if ~isrow(Tr)
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(n18), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(n18) && isrow(n18)
            n18 = n18';
            n0 = n0';
        elseif ~isvector(n18) && size(n18,2) ~= nTheta
            error('If n18 is a matrix, size(n18,2) has to be size(Tr,2)')
        end
    end
    
    if ~isfield(p, 'J')
        p.J = size(n18,2);                                                  % number of considered reactions (-)
    end
    
    A = p.p_est(1);
    k = p.p_est(2);
    n1 = p.p_est(3);
    if numel(p.p_est) > 3
        n2 = p.p_est(4); 
    else
        n2 = 0;
    end

    if nargin < 7
        d = const.d;
    end
    
    r = const.r;                                                            % radius of support material (m)
    rPe = r + d;                                                            % radius of pellet with active layer (m)
  
    HG = n18./(n18 + n0);                                                   % hydrogenation grade (-)    
    m = n18*const.M18 + n0*const.M0;                                        % mass of liquid phase (g)
    rho = getRho(const, pr, Tr, HG);                                        % density of liquid mixture, wo H2 (g l^-1)
    V = m./rho;                                                             % volume of liquid phase (l)
    
    if isfield(p, 'LRxz')
        epsBed = const.epsBed;
    else
        epsBed = const.epsBed;
    end
    VPe = p.Vc.*(1-epsBed)/1000;                                            % volume of pellets (l)
    Vbed = V + VPe;                                                         % volume of bed (l)

    Trb = Tr;                                                               % reaction temperature in bulk (K) 
    c18b = n18./V;                                                          % molar concentration of H18-DBT in bulk (mol l^-1)    
    switch p.rate
        case 0
            frate = const.M18;                                              % conversion factor from c18 to partial density rho18 (g mol^-1)                             
        case 1
            frate = 1;                                                      % conversion factor from c18 to molar concentration c18 in mol l^-1 (-)
        case 2        
%             frate = V./(n18 + n0);                                          % conversion factor from c18 to molar fraction x18 (-)                
            frate = const.M18./rho;                                         % conversion factor from c18 to mass fraction w18 (-)
    end 

    fkr = 1;                                                                % conversion factor for kr
%     fkr = p.nL0.^-1;
%     fkr = p.xPt;
%     fkr = p.xPt./(1-mean(p.epsBed))*100;
%     fkr = p.mc./p.V0*1000;                                                  % factor from AS & PP
%     fkr = p.mc./p.m0;
%     fkr = p.Vc./p.V0*1000;
%     fkr = p.V0./(p.V0+VPe*1000);
%     fkr = p.mc./(p.V0+p.Vc.*(1-epsBed))*1000;
%     fkr = p.mc./(p.nL0.*(HG.*const.M18+(1-HG).*const.M0)+p.mc);
%     if isfield(p, 'mPt')
%         fkr = p.mPt./p.V0*1000;                                             
%     else
%         mPt = p.xPt.*p.nL0.*const.MPt;
%         fkr = mPt./p.V0*1000; 
%     end

    fkr2 = 1;                                                               % conversion factor for n
%     fkr2 = (p.xPt.*p.nL0.*10000);           
%     fkr2 = (p.nL0.*const.M18./getRho(const, pr, Tr, 1)+VPe)./p.Vr*1000;
%     fkr2 = p.Vc;
%     fkr2 = p.xPt.*(p.V0+p.Vc.*(1-epsBed))*100;
%     fkr2 = p.V0./(p.V0+p.Vc.*(1-epsBed));
    
    rateB4 = myRate(Trb, c18b.*frate);                                      % reaction rate @ bulk conditions (mol l_bed^-1 s^-1)
    rateB3 = rateB4.*Vbed./VPe;                                             % reaction rate @ bulk conditions (mol l_c^-1 s^-1)
    reff4 = rateB4;                                                         % effective reaction rate (mol l_bed^-1 s^-1)
    etac = ones(size(rateB4));                                              % effectivenes factor of catalyst (-)
    P = [];                                                                 % concentration and temperature profile inside catalyt (mol l^-1) and (K)
    if p.ir == 1 && any(n0 > 0)                                               
        if (rPe-r)*1e6 < 55
            dr = 5;
        elseif (rPe-r)*1e6 < 105
            dr = 10;
        elseif (rPe-r)*1e6 < 155
            dr = 10;
        else
            dr = 20;
        end
%         N = round((rPe-r)*1e6/dr);                                          % number of volume elements (-)
%         rsol = linspace(r, rPe, N+1);
        rsol = [r:dr*1e-06:r+(round((rPe-r)*1e06,-1)-5)*1e-06 rPe];
        rinit = (linspace(0, rPe, 10)./(rPe)).^(1/3);                        % cubic root spacing (dimensionsless)
        rinit = rinit.*(rPe);

        dHr = getDHr(const, pr, Trb);                                       % reaction enthalpy for dehydrogenation of H18-DBT (kJ mol^-1)    
        DPe18 = getDPe(const, p, Trb, HG);                                  % binary diffusion coefficient of H18-DBT in liquid mixture, wo H2 (m^2 s^-1)
        [~, lambdaPed] = getLambdaPe(const, p, pr, Tr, HG);                 % effective heat conductivity of active layer (catalyst) (W m^-1 K^-1)
        
        EA = A*(const.R.*p.Tm(1));
        p.AndersonH = 3*const.R*Trb.^2.*lambdaPed./ ...                     % Anderson criterion with rate in (mol l_c^-1 s^-1), for >= 1 no heat transport limitations in pellet (-)
            ((1000*rateB3).*(1000*dHr).*EA.*(2*rPe).^2);   
        if any(p.AndersonH < 1) && p.irT == 0 && ~isfield(p, 'warnA')
            p.warnA = fprintf(['Internal HT limitation should ' ...
                'be considered (set p.irT = 1)!\n']);
        end
        p.MearsH = (p.BiH.*lambdaPed.*Trb)./ ...                            % Mears criterion, for Mears > 1 no heat transport limitation in film  
            ((1000*dHr).*(1000*rateB3).*rPe.^2).*0.15./(EA./(const.R.*Trb));    
        p.MearsM = (p.BiM.*DPe18.*c18b)./(rateB3.*rPe.^2).* ...             % Mears criterion, for MearsM > 1 no mass transport limitation in film 
            0.15./abs(n1 + n2.*pr);
        if any(p.MearsH <= 1) 
            if (p.film == 0 || p.BiH >= 1e03) && ~isfield(p, 'warnMH')
               p.warnMH = fprintf(['HT limitation in film should ' ...
                   'be considered (set p.film = 1 & BiH << 1000)!\n']);
            end
        elseif any(p.MearsM <= 1) 
            if (p.film == 0 || p.BiM >= 1e03) && ~isfield(p, 'warnMM')
                p.warnMM = fprintf(['MT limitation in film should ' ...
                'be considered (set p.film = 1 & BiM << 1000)!\n']);
            end
        end

        if p.irT == 0
            solinit = bvpinit(rinit, ...
                [c18b.*frate zeros(1,p.J)]);
        else
            solinit = bvpinit(rinit, ...
                [c18b.*frate zeros(1,p.J) Trb zeros(1,p.J)]);
        end

        S = diag(repmat([zeros(1,p.J) -2.*ones(1,p.J)], 1, ...              % treatment of singularity for pellet (sphere)
            numel(solinit.yinit)/(2*p.J)));
        options = bvpset('Vectorized', 'on', 'SingularTerm', S);

        sol = bvp4c(@ODE, @BC, solinit, options);
        Y = deval(sol, rsol);

        c18 = Y(1:p.J,:)'./frate;                                           % molar concentration of H18-DBT (mol l^-1)
        c18(rsol < r) = 0;
        if p.irT == 0
            Tr = Trb; 
            Trs = Trb;
        else
            Tr = Y(2*p.J+1:3*p.J,:)';
            Trs = Tr(end,:);                                                % temperature @ catalyst surface (K)
        end
        
        reff4 = 1./(4/3*pi*((rsol(end)).^3-(rsol(1)).^3)).* ...             % effective reaction rate (mol l_bed^-1 s^-1)
            trapz((rsol), myRate(Tr, c18.*frate)*4*pi.*rsol'.^2, 1);
        reff3 = reff4.*Vbed./VPe;                                           % effective reaction rate (mol l_c^-1 s^-1)
        if size(Tr,1) ~= size(rsol',1)
            Tr = repmat(Tr,size(rsol'));
        end

        etac = reff3./rateB3;                                               % effectivenes factor of catalyst (-)
        P = [rsol' c18 Tr];                                                 % concentration and temperature profile inside catalyt (mol l^-1) and (K)

        dTmax = getDTmax(const, p, pr, Trs, HG);                            % maximum possible temperature difference inside the active layer (K)
        if any(Trs - Tr(1,:) > dTmax) && ~isfield(p, 'warnDT')
            p.warnDT = fprintf(['Maximum possible temperature ' ...
                'difference inside pellet exceeded (vary p.lambdaPe)!\n']);
        end
    end
   

    %% nested functions
    function rate = myRate(Tr, conc18)                                      % conc18 either rho18 (p.rate=0) or c18 (p.rate=1)
        n = n1 + n2.*pr;                                                    % reaction order (-)
        kr = k.*fkr.*exp(A.*(1-p.Tm./Tr));                                  % frequency of collision (mol^(1-n) s^-1  l^(n-1) frate^-n)
        rate = kr.*conc18.^(n.*fkr2);                                       % reaction rate (mol l_bed^-1 s^-1)

%         concH2 = pr*1e05./(const.R.*Tr).*frate;                             % cocentration of H2 in gas phase (ideal gas) (mol s^-1)
%         kr = k.*fkr.*exp(A.*(1-p.Tm./Tr));                                  % frequency of collision (mol^(1-n) s^-1  l^(n-1) frate^-n)
%         rate = kr.*conc18.^(n1.*fkr2).*concH2.^n2;                          % reaction rate (mol l_bed^-1 s^-1)        
    end

    function dYdr = ODE(R, Y)
        if p.irT == 0
            TR = Trb;
        else
            TR = Y(2*p.J+1:3*p.J,:)';
        end
        
        rate3 = myRate(TR, Y(1:p.J,:)').*Vbed./VPe;                         % reaction rate inside pellet (mol l_c^-1 s^-1)
        rate3(R <= r) = 0;

        dYdr = [Y(p.J+1:2*p.J,:); ...
                -(const.nu(1).*frate.*rate3./DPe18)'];               
        if ~(p.irT == 0)
            dYdr = [dYdr; ...
                    Y(3*p.J+1:4*p.J,:); ...
                    ((1000.*dHr).*(1000.*rate3)./lambdaPed)'];
        end
    end

    function res = BC(Y0, YPe)                                              % Y0 is boundary at r=0 (center) and YPe is boundary at r = rPe (surface of pellet)
        if p.film == 0
            res = [YPe(1:p.J,:) - (c18b.*frate)'; ...
                    Y0(p.J+1:2*p.J,:)];
            if ~(p.irT == 0) 
                res = [res; ...
                        YPe(2*p.J+1:3*p.J,:) - Trb'; ...
                        Y0(3*p.J+1:4*p.J,:)];
            end   
        else
            res = [YPe(p.J+1:2*p.J,:) - ...
                        p.BiM./rPe.*((c18b.*frate)' - YPe(1:p.J,:)); ...
                    Y0(p.J+1:2*p.J,:)];
            if ~(p.irT == 0) 
                res = [res; ...
                        YPe(3*p.J+1:4*p.J,:) - ...
                            p.BiH./rPe.*(Trb' - YPe(2*p.J+1:3*p.J,:)); ...
                        Y0(3*p.J+1:4*p.J,:)];
            end   
        end

    end
end