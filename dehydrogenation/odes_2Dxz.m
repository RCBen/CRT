function dXdx = odes_2Dxz(x, X, const, p)
    %
    % function for calculating the derivatives
    
    % disassemble the vector of states
    rhou0 = X(0*p.nz+1:1*p.nz)';                                            % mass flux in x-direction (kg m^-2 s^-1)
    w18 = X(1*p.nz+1:2*p.nz)';                                              % mass fraction of H18-DBT (-)
    Tr = X(2*p.nz+1:3*p.nz)';                                               % reaction temperature (K)
    Tw = p.Tw(dsearchn(p.x', x));                                           % wall temperature (K)
    
    Mn = (w18./const.M18 + (1-w18)./const.M0).^-1;                          % mean molar mass of liquid mixture (g mol^-1)
    HG = w18./const.M18.*Mn;                                                % hydrogenation grade (-)
    
    % substance properties    
    rho = getRho(const, [], Tr, HG);                                        % density of liquid mixture, wo H2 (g l^-1)
    rhoH2Id = ((p.pr*1e+05).*const.MH2./(const.R.*Tr))*1e-03;               % density of ideal gas H2 (g l^-1)

    [DEff18, DEff0] = getDEff(const, p, Tr, HG);                            % binary diffusion coefficient of H18-DBT & H0-DBT in liquid mixture, wo H2 (m^2 s^-1)
    [Qm, sigma18, sigma0, sigmaH2] = myQm(const, p, p.pr, Tr, HG);          % mass source & species rate for consumption of H18-DBT & production of H0-DBT & H2 (g l_bed^-1 s^-1)
    
    if p.iso == 0
        cp = (getCp(const, [], Tr, HG)./Mn)*1e06;                           % isobaric heat capacity of liquid mixture, wo H2 (J kg^-1 K^-1)
        lambdaEff = getLambdaEff(const, p, [], Tr, HG);                     % effective heat conductivity of catalyst bed (W m^-1 K^-1)   
        [h18, h0, hH2] = getH(const, p.pr, Tr);                             % molar enthalpy of H18-DBT & H0-DBT & H2 (kJ mol^-1)
        h18m = (h18./const.M18)*1e06;                                       % enthalpy of H18-DBT (J kg^-1)
        h0m = (h0./const.M0)*1e06;                                          % enthalpy of H0-DBT (J kg^-1)
        hH2m = (hH2./const.MH2)*1e06;                                       % enthalpy of H2 (J kg^-1)       
    end
    
    u0 = rhou0./rho;                                                        % superficial velocity of liquid in x-direction (m s^-1)
    u0(p.z == 0) = 0;                                                       % superficial velocity of liquid in x-direction @ wall (m s^-1)

    w0H2 = cumtrapz(p.z, -Qm)./rhoH2Id./p.epsBed;                           % superficial velocity of H2 in z-direction (m s^-1)
    
    if p.iso == 0
        kamb = p.kamb(Tr);                                                  % ambient heat transfer coefficient (W m^-2 K^-1)    
                                                                            % kamb is referred only to reaction volume (without reactor and heating plate)
%         kamb = zeros(size(Tr));                                                                            
        kw = myKw(const, p);                                                % overall heat transfer coefficient @ wall (W m^-2 K^-1)
    end

    % calculating the derivatives (central difference)
    rhow18 = rho.*w18;
    rhow0 = rho.*(1-w18);
    i = 1;
    drhou0dx(i) = Qm(i);
    dw18dx(i) = 1./rhou0(i).*(sigma18(i) - w18(i).*Qm(i) + ...
        ((DEff18(i+1) + DEff18(i))/2.*(rhow18(i+1)-rhow18(i)) - ...
        0)./p.dz^2);
    if p.iso == 0
        dTrdx(i) = 1./(rhou0(i).*cp(i)).* ...
            (-(sigma18(i).*h18m(i) + sigma0(i).*h0m(i) + sigmaH2(i).*hH2m(i)) + ...
            ((lambdaEff(i+1) + lambdaEff(i))/2.*(Tr(i+1)-Tr(i)) - ...
            kw.*(Tr(i)-Tw).*p.dz)./p.dz^2 + ...
            (DEff18(i).*((rhow18(i+1) + rhow18(i))/2.* ...
            (const.g.*p.dz + (h18m(i+1)-h18m(i))) - ...
            (rhow18(i) + rhow18(i))/2.* ...
            (const.g.*p.dz + 0)) + ...
            DEff0(i).*((rhow0(i+1) + rhow0(i))/2.* ...
            (const.g.*p.dz + (h0m(i+1)-h0m(i))) - ...
            (rhow0(i) + rhow0(i))/2.* ...
            (const.g.*p.dz + 0)))./p.dz.^2 - ...
            kamb(i).*p.SvRxz.*(Tr(i) - p.Tamb));
    end
    
    i = 2:p.nz-1;
    drhou0dx(i) = Qm(i);
    dw18dx(i) = 1./rhou0(i).*(sigma18(i) - w18(i).*Qm(i) + ...
        ((DEff18(i+1) + DEff18(i))/2.*(rhow18(i+1)-rhow18(i)) - ...
        (DEff18(i) + DEff18(i-1))/2.*(rhow18(i)-rhow18(i-1)))./p.dz^2);
    if p.iso == 0
        dTrdx(i) = 1./(rhou0(i).*cp(i)).* ...
            (-(sigma18(i).*h18m(i) + sigma0(i).*h0m(i) + sigmaH2(i).*hH2m(i)) + ...
            ((lambdaEff(i+1) + lambdaEff(i))/2.*(Tr(i+1)-Tr(i)) - ...
            (lambdaEff(i) + lambdaEff(i-1))/2.*(Tr(i)-Tr(i-1)))./p.dz^2 + ...
            (DEff18(i).*((rhow18(i+1) + rhow18(i))/2.* ...
            (const.g.*p.dz + (h18m(i+1)-h18m(i))) - ...
            (rhow18(i) + rhow18(i-1))/2.* ...
            (const.g.*p.dz + (h18m(i)-h18m(i-1)))) + ...
            DEff0(i).*((rhow0(i+1) + rhow0(i))/2.* ...
            (const.g.*p.dz + (h0m(i+1)-h0m(i))) - ...
            (rhow0(i) + rhow0(i-1))/2.* ...
            (const.g.*p.dz + (h0m(i)-h0m(i-1)))))./p.dz.^2 - ...
            kamb(i).*p.SvRxz.*(Tr(i) - p.Tamb));
    end

    i = p.nz;
    drhou0dx(i) = Qm(i);
    dw18dx(i) = 1./rhou0(i).*(sigma18(i) - w18(i).*Qm(i) + ...
        (0 - ...
        (DEff18(i) + DEff18(i-1))/2.*(rhow18(i)-rhow18(i-1)))./p.dz^2);
    if p.iso == 0
        dTrdx(i) = 1./(rhou0(i).*cp(i)).* ...
            (-(sigma18(i).*h18m(i) + sigma0(i).*h0m(i) + sigmaH2(i).*hH2m(i)) + ...
            (0 - ...
            (lambdaEff(i) + lambdaEff(i-1))/2.*(Tr(i)-Tr(i-1)))./p.dz^2 + ...
            (DEff18(i).*((rhow18(i) + rhow18(i))/2.* ...
            (const.g.*p.dz + 0) - ...
            (rhow18(i) + rhow18(i-1))/2.* ...
            (const.g.*p.dz + (h18m(i)-h18m(i-1)))) + ...
            DEff0(i).*((rhow0(i) + rhow0(i))/2.* ...
            (const.g.*p.dz + 0) - ...
            (rhow0(i) + rhow0(i-1))/2.* ...
            (const.g.*p.dz + (h0m(i)-h0m(i-1)))))./p.dz.^2 - ...
            kamb(i).*p.SvRxz.*(Tr(i) - p.Tamb));
    end
    
    if p.iso == 1
        dTrdx = zeros(size(p.z));
    end
    dXdx = [drhou0dx dw18dx dTrdx]';

    %% nested functions
    function kw = myKw(const, p)
        % see DOI: 10.1007/978-3-642-19981-3 (pp. 27 & 805 ff & 1521 & 1527)
        %    ISBN: 3-527-31000-2 (pp. 51 & 242)
        
        Tm = mean([Tr(p.z == 0), Tw]);
        x18m = HG(p.z == 0);
        
        rhom = getRho(const, [], Tm, x18m);                                 % density of liquid mixture, wo H2 (g l^-1)
        etam = getEta(const, [], Tm, x18m)*1e-03;                           % dynamic viscosity of liquid mixture, wo H2 (Pas)
        cpm = getCp(const, [], Tm, x18m)*1e06./ ...                         % isobaric heat capacity of liquid mixture, wo H2 (J kg^-1 K^-1)
            (x18m.*const.M18+(1-x18m).*const.M0);   
        lambdam = getLambda(const, [], Tm, x18m);                           % heat conductivity of mixture, wo H2 & ideal gas H2 (W m^-1 K^-1)  
             
        Pr = etam.*cpm./lambdam;                                            % Prandtl number (-)
        switch p.hw
            case 1
                %1 - plate (Tw = const.)
                Uext = mean(u0)./p.epsBed;                                  % exterior overflow velocity (far from boundary layer) of liquid in x-direction (m s^-1)
                               
                Rex = rhom.*Uext.*x./etam;                                  % local Reynolds number (-)
                if all(Rex(:) <= 1e07)
                    Nulamx = 0.332.*sqrt(Rex).*Pr.^(1/3);                   % local Nusselt number for laminar flow (-)
                    Nuturbx = 0.0296.*Rex.^0.8.*Pr./ ...                    % local Nusselt number for turbulent flow (-)
                        (1+2.185.*Rex.^(-0.1).*(Pr-1));                      
                    Nux = sqrt(Nulamx.^2+Nuturbx.^2);                       % local Nusselt number (-)
                else
                    fprintf(['Rex > 1e07, thus calculation' ...
                'not valid (set p.hw to 2, 3 or 4)!\n']);
                end
                
                hw = Nux.*lambdam./max(x, sqrt(eps));
            case 2
                %2 - packed bed (Tw = const.)
                U0 = mean(u0);                                              % mean superficial velocity of liquid in x-direction (m s^-1)
                
                Rep0 = rhom.*U0.*(2*const.rd)./etam;                        % particle Reynolds number (-)
                Pe0 = Rep0.*Pr;                                             % Peclet number (-)
                if all(Pe0(:) > 1) && all(Pe0(:) < 1e04)
                    lambdaEffm = getLambdaEff(const, p, [], Tm, x18m);      % effective heat conductivity of catalyst bed (W m^-1 K^-1)
                    Nu = (1.3 + 5./(p.DhRxz/(2*const.rd))).* ...            % Nusselt number (-)
                        lambdaEffm./lambdam + 0.19*Rep0.^0.75.*Pr.^(1/3);
                end

                hw = Nu.*lambdam./(2*const.rd);
            case 3
                %3 - bubble column (Tw = const.)
                w0G = mean(w0H2);                                           % mean superficial velocity of gas in z-direction (m s^-1)
                
                Reb0G = rhom.*w0G*p.db./etam;                               % bubble Reynolds number (-)
                FrG = w0G.^2./(const.g*p.db);                               % Froude number of gas (-)
                St = 0.1*(Reb0G.*FrG.*Pr.^2).^-0.25;                        % Stanton number (-)
                
                hw = St.*(rhom.*cpm.*w0G);
            case 4
                %4 - Inf (T(z=0) = Tw)
                hw = 1e12;
        end

        delta = 3.5e-03;                                                    % distance between wall and thermocouple (m)        
        lambdaSS = 16.4658 - 0.0050*Tm;                                     % heat conductivity of stainless steel (316Ti, 1.4571)
        kw = (1./hw + delta./lambdaSS).^-1;                                 % overall heat transfer coefficient @ wall (W m^-2 K^-1);
    end
end