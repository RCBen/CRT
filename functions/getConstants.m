function const = getConstants
    % const = getConstants
    %
    % physical constants & 
    % constant properties of reaction, catalyst (solid), H2 (gas) & LOHC (liquid)
    % 
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 141 ff)
    %          10.1021/acs.iecr.5b01840
    %          10.1016/j.fluid.2017.01.017
    %    ISBN: 0-07-0011682-2 (pp. 2.23 & 4.8 & 7.7 f & A.19)
    %     LOT: CM23001 (Pt/Alox, Hydrogenious)
    %     URL: http://webbook.nist.gov/chemistry/
    
    %% physical constant
    const.g = 9.8066;                                                       % acceleration of gravity (m s^-2)
    const.R = 8.3143;                                                       % gas constant (J mol^-1 K^-1)
    const.NA = 6.0221e23;                                                   % Avogadro constant (mol^-1)
    const.kB = const.R/const.NA;                                            % Boltzmann constant (J K^-1)
    const.sigma = 5.67*1e-08;                                               % Stefan-Boltzmann constant (W m^-2 K^-4)
    
    %% reaction
    % see DOI: 10.1021/acs.iecr.5b01840
    %          10.1016/j.fluid.2017.01.017
    %     URL: http://webbook.nist.gov/chemistry/
    const.nu = [-1; 1; 9];                                                  % stoichiometric coefficients of H18-DBT & H0-DBT & H2 (-)
    const.p0 = 1.01325;                                                     % pressure @ standard state (bar)
    const.T0 = 298.15;                                                      % temperature @ standard state (K)
    const.dHb18 = -475.6;                                                   % standard enthalpy of formation for H18-DBT (kJ mol^-1)
    const.dHb0 = 112.9;                                                     % standard enthalpy of formation for H0-DBT (kJ mol^-1)
    const.dHbH2 = 0;                                                        % standard enthalpy of formation for H2 (kJ mol^-1)
    const.dHr0 = sum([const.dHb18; const.dHb0; const.dHbH2].*const.nu) + ...% standard reaction enthalpy for dehydrogenation of H18-DBT (wo uncertainty) (kJ mol^-1)
        0*sqrt(5.6^2 + 7.6^2);
    const.dSr0 = (3*(220.96 + 3*130.68 - 247.90) + 9*7.8)/1000;             % standard reaction entropy for dehydrogenation of H18-DBT (kJ mol^-1 K^-1)
                                                                            % 3x standard reaction entropy for dehydrogenation of MCH (w uncertainty)
    const.dGr0 = (const.dHr0 - const.T0*const.dSr0);                        % standard gibbs free energy for dehydrogenation of H18-DBT (kJ mol^-1)
    const.lnK0 = -(1000*const.dGr0)/(const.R*const.T0);                     % logarithm of standard equilibrium constant (-)
    
    %% catalyst (solid)
    % see LOT: CM23001 (Pt/Alox, Hydrogenious)
    %     URL: http://webbook.nist.gov/chemistry/
    const.MPt = 195.084;                                                    % molar mass of Pt (g mol^-1)
    const.wPtc = 0.005;                                                     % mass fraction of Pt regarding to mass of cat (-)
    const.lambdaPt = 73;                                                    % heat conductivity of Pt (W m^-1 K^-1)   
    
    const.d = 95e-6;                                                        % thickness of active layer (m)    
    const.r = 0.0015;                                                       % radius of pellet (m)       
%     const.r = 0.0015 - const.d;                                             % radius of pellet (m)  
    const.rd = const.r + const.d;                                           % radius of pellet with active layer (m)
    const.V = 4/3*pi*(const.rd^3)*1000;                                     % volume of pellet (l)
    const.phiPt = 1-(const.r/const.rd).^3;                                  % volume fraction of catalyst with Pt
    
    const.dpore = 3.161*1e-09;                                              % average pore diameter (m)
%     const.dpore = mean([5.564, 6.155])*1e-09;                               % average pore diameter (m)
    
%     const.rhoPe = round(0.9953/(37*4/3*pi*(const.rd)^3)*1e-03);             % solid density of porous catalyst (g l^-1)
    const.rhoPe = 1486;                                                     % solid density of porous catalyst (Rohdichte) (g l^-1)
    const.rhoAlox = 3900;                                                   % solid density of alumina (Reindichte) (g l^-1)
    const.epsPe = 1-const.rhoPe/const.rhoAlox;                              % internal porosity of catalyst (-)
    
    const.rhoBed = 800;                                                     % bed density of porous domain (Schüttdichte) (g l^-1)
    const.epsBedInf = 0.40;                                                 % unbounded bed porosity(-)
    const.epsBed = 1-const.rhoBed/const.rhoPe;                              % mean bed porosity (-)
    const.rhoBedInf = (1-const.epsBedInf)*const.rhoPe;                      % mean bed density (g l^-1)
    
    %% H2 (gas, classical (wH2 = 0) vs. quantum (wH2 < 0))
    % see DOI: 10.1002/aic.690190120
    %    ISBN: 0-07-0011682-2 (p. 4.8 & A.19)    
    %     URL: http://webbook.nist.gov/chemistry/
    const.MH2 = 2*1.00794;                                                  % molar mass of H2 (g mol^-1)
    const.muH2 = 0;                                                         % dipole moment of H2 (D)
    
    const.TcH2 = 43.60;                                                     % critical temperature of H2 (classical) (K)
    const.pcH2 = 20.5;                                                      % critical pressure of H2 (classical) (bar)
    const.vcH2 = 0.0515;                                                    % critical molar volume of H2 (classical) (l mol^-1)
    const.ZcH2 = (const.pcH2*1e+05)*(const.vcH2/1000)/ ...                  % critical compressibility factor of H2 (classical) (-)
        (const.R*const.TcH2);   
    const.wH2 = 0;                                                          % acentric factor of H2 (classical) (-)
    
    const.TcH2q = 32.98;                                                    % critical temperature of H2 (quantum) (K)
    const.pcH2q = 12.93;                                                    % critical pressure of H2 (quantum) (bar)
    const.vcH2q = 0.0642;                                                   % critical molar volume of H2 (quantum) (l mol^-1)
    const.ZcH2q = (const.pcH2q*1e+05)*(const.vcH2q/1000)/ ...               % critical compressibility factor of H2 (quantum) (-)
        (const.R*const.TcH2q);
    const.wH2q = -0.2170;                                                   % acentric factor of H2 (quantum) (-)    
    
    const.pn = 1.01325;                                                     % pressure @ normal state (bar)
    const.Tn = 293.15;                                                      % temperature @ normal state (K)
    const.rhoH2n = 0.083755;                                                % density of H2 @ normal state (g l^-1)
    const.LHVH2 = 120*1000;                                                 % lower heating value of H2 (kJ kg^-1)
    
    %% LOHC (liquid)
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 141 ff)
    %    ISBN: 0-07-0011682-2 (pp. 2.23 & 7.7 f)
    %     URL: http://webbook.nist.gov/chemistry/
    const.M18 = 21*12.0107 + 38*1.00794;                                    % molar mass of H18-DBT (g mol^-1)
    const.mu18 = 0;                                                         % dipole moment of H18-DBT (D)
    const.dmol18 = 1.3*1e-09;                                               % molecule diameter (m)
    
    const.Tc18 = 181.128*log(15*3.492 + 1*1.6781 + 5*4.033 + 3*0.8479);     % critical temperature of H18-DBT (K)
    const.pc18 = (0.10022 + ...                                             % critical pressure of H18-DBT (bar)
        (15*0.010558 + 1*0.019904 + 5*0.001315 + 3*0.002257))^-2 + ...
        1.3705;                                                   
    const.vc18 = (15*0.05576 + 1*0.07504 + 5*0.03153 + 3*0.01636) - ...     % critical molar volume of H18-DBT (l mol^-1)
        0.00435;                                                    
    const.Zc18 = (const.pc18*1e+05)*(const.vc18/1000)/ ...                  % critical compressibility factor of H18-DBT (-)
        (const.R*const.Tc18);
    
    const.Tb18 = 621.6532;                                                  % normal boiling point of H18-DBT (K)
    Tbr18 = const.Tb18/const.Tc18;                                          % reduced temperature of H18-DBT (-)
    tau18 = 1 - Tbr18;
    f180 = (-5.97616*tau18 + 1.29874*tau18^1.5 - ...
        0.60394*tau18^2.5 - 1.06841*tau18^5)/Tbr18;
    f181 = (-5.03365*tau18 + 1.11505*tau18^1.5 - ...
        5.41217*tau18^2.5 - 7.46628*tau18^5)/Tbr18;
    const.w18 = -(log(const.pc18/1.01325) + f180)/f181;                     % acentric factor of H18-DBT (-)

    const.M0 = 21*12.0107 + 20*1.00794;                                     % molar mass of H0-DBT (g mol^-1)
    const.mu0 = 0;                                                          % dipole moment of H0-DBT (D)
    
    const.Tc0 = 181.128*log(2*14.6409 + 13*3.7337 + 2*10.3239 + 1*8.2130);  % critical temperature of H0-DBT (K)
    const.pc0 = (0.10022 + ...                                              % critical pressure of H0-DBT (bar)
        (2*0.002136 + 13*0.007542 + 2*0.012200 + 1*0.019360))^-2 + ...
        1.3705;
    const.vc0 = (2*0.03985 + 13*0.04215 + 2*0.10099 + 1*0.10364) - ...      % critical molar volume of H0-DBT (l mol^-1)
        0.00435;
    const.Zc0 = (const.pc0*1e+05)*(const.vc0/1000)/ ...                     % critical compressibility factor of H0-DBT (-) 
        (const.R*const.Tc0);         
    
    const.Tb0 = 640.3431;                                                   % normal boiling point of H0-DBT (K)
    Tbr0 = const.Tb0/const.Tc0;                                             % reduced normal boiling point of H0-DBT (-)
    tau0 = 1 - Tbr0;
    f00 = (-5.97616*tau0 + 1.29874*tau0^1.5 - ...
        0.60394*tau0^2.5 - 1.06841*tau0^5)/Tbr0;
    f01 = (-5.03365*tau0 + 1.11505*tau0^1.5 - ...
        5.41217*tau0^2.5 - 7.46628*tau0^5)/Tbr0;
    const.w0 = -(log(const.pc0/1.01325) + f00)/f01;                         % acentric factor of H0-DBT (-)
    
    vB = [const.M18; const.M0]./...
        getRho(const, [], [const.Tb18, const.Tb0], [1, 0])';
    const.vB18 = vB(1);                                                     % molar volume of H18-DBT @ Tb (l mol^-1)
    const.vB0 = vB(2);                                                      % molar volume of H0-DBT @ Tb (l mol^-1)
     
%     % see DOI: 10.1007/978-3-642-19981-3 (pp. 149 f)
%     %    ISBN: 0-07-0011682-2 (pp. 4.38 f)
%     VR180 = 1 - 1.52816.*(1-(const.Tb18./const.Tc18)).^(1/3) + ...
%         1.43907.*(1-(const.Tb18./const.Tc18)).^(2/3) - ...
%         0.81446.*(1-(const.Tb18./const.Tc18)) + ...
%         0.190454.*(1-(const.Tb18./const.Tc18)).^(4/3);
%     VR18delta = (-0.296123 + 0.386914.*(const.Tb18./const.Tc18) - ...
%         0.0427458.*(const.Tb18./const.Tc18).^2 - ...
%         0.0480645.*(const.Tb18./const.Tc18).^3)./...
%         ((const.Tb18./const.Tc18) - 1.00001);
%     const.vB18 = const.vc18.*VR180*(1-const.w18.*VR18delta);              % molar volume of H18-DBT @ Tb (l mol^-1)

%     VR00 = 1 - 1.52816.*(1-(const.Tb0./const.Tc0)).^(1/3) + ...
%         1.43907.*(1-(const.Tb0./const.Tc0)).^(2/3) - ...
%         0.81446.*(1-(const.Tb0./const.Tc0)) + ...
%         0.190454.*(1-(const.Tb0./const.Tc0)).^(4/3);
%     VR0delta = (-0.296123 + 0.386914.*(const.Tb0./const.Tc0) - ...
%         0.0427458.*(const.Tb0./const.Tc0)^2 - ...
%         0.0480645.*(const.Tb0./const.Tc0)^3)/...
%         ((const.Tb0./const.Tc0) - 1.00001);
%     const.vB0 = const.vc0.*VR00*(1-const.w0.*VR0delta);                   % molar volume of H0-DBT @ Tb (l mol^-1) 

    const.wH218 = const.nu(3)*const.MH2/const.M18;                          % mass fraction of releasable H2 in H18-DBT (-)
end