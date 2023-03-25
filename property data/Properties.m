% function PropertiesFit(prop, hg, mix, excl)
%
% property prop             2   nD
%                           3   cp
%                           4   rho
%                           5   eta
%                           6   pLV
%                           7   hLV
%                           8   lambda
%                           9   omega
%
% hydrogenation grade hg    1   NMR
%                           2   EA
%
% mixture mix               1   reaction
%                           2   artificial
%
% exclude excl              0   no exclusion
%                           1   exclude larg
%
    hg = 1;
%   mix = 1; prop = 2; excl = 0;
%    mix = 2; prop = 2; excl = 0;
%    mix = 1; prop = 3; excl = 1;
%    mix = 1; prop = 4; excl = 0;
%    mix = 2; prop = 4; excl = 0;
%    mix = 1; prop = 5; excl = 0;
%    mix = 2; prop = 5; excl = 0;
    mix = 1; prop = 6; excl = 1;
%     mix = 2; prop = 6; excl = 1;
    
    method = {'NMR', 'EA'};
    mixture = {'r', 'm'};
    sheet = {'HG', 'nD', 'cp', 'rho', 'eta', 'pLV', 'hLV', 'lambda', 'omega'};
    unit = {'-', '-', 'J g^{-1} K^{-1}', 'g cm^{-3}', 'mPa s', 'mbar', 'kJ mol^{-1}', 'W m^{-1} K^{-1}', 'mN m^{-1}'};
    
    %% hydrogenation grade
    HGmat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{1});
    sampler = HGmat(2:17,1);
    samplem = [1,3,4,6,7,8,16,11,12,14]';
    nrm = [1,3,4,6,7,8,11,13,14,16];
    if hg == 1
        HGr = HGmat(2:17,2);
        HGm = HGmat(nrm+1,7);
    elseif hg == 2
        HGr = HGmat(2:17,3);
        HGm = HGmat(nrm+1,8);
    end
    HGr(1) = 0;
    HGr(end) = 1;
    HGm(1) = 0;
    HGm(end) = 1;
    
    %% property
    if prop == 2
        nDmat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{prop});
        Theta = nDmat(1,2:7);                                               % temperature / °C
        if mix == 1
            nD = nDmat(2:17,2:7);                                           % refractive index of reaction mixtures / nD
        elseif mix == 2
            nD = nDmat(nrm+1, 11:16);                                       % refractive index of artificial mixtures / nD
        end
    elseif prop == 3
        cpmat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{prop});
        Theta = cpmat(1,2:9);                                               % temperature / °C
        if mix == 1
            cp = cpmat(2:17,2:9);                                           % heat capacity of reaction mixtures / J g^-1 K^-1
%         elseif mix == 2
%             cp = cpmat(nrm+1, 13:20);                                     % heat capacity of artificial mixtures / J g^-1 K^-1
        end
    elseif prop == 4
        rhomat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{prop});
        Theta = rhomat(1,2:7);                                              % temperature / °C
        if mix == 1
            rho = rhomat(2:17,2:7);                                         % density of reaction mixtures / g cm^-3
        elseif mix == 2
            rho = rhomat(nrm+1, 11:16);                                     % density of artificial mixtures / g cm^-3
        end
    elseif prop == 5
        etamat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{prop});
        Theta = etamat(1,4:20);                                             % temperature / °C
%         Theta = etamat(1,17:20);
        if mix == 1
            eta = etamat(2:17,4:20);                                        % viscosity of reaction mixtures / mPa s
%             eta = etamat(2:17, 17:20);
        elseif mix == 2
            eta = etamat(nrm+1, 26:42);                                     % viscosity of artificial mixtures / mPa s
%             eta = etamat(nrm+1, 39:42);
        end
    elseif prop == 6
        pLVmat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{prop});
        if mix == 1
            Theta = pLVmat(2:17,2:15);                                      % temperature / °C
            pLV = pLVmat(23:38,2:15);                                       % viscosity of reaction mixtures / mPa s
        elseif mix == 2
            Theta = pLVmat(nrm+1,20:33);                                    % temperature / °C
            pLV = pLVmat(nrm+22, 20:33);                                    % viscosity of artificial mixtures / mPa s
        end
    elseif prop == 9
        omegamat = readmatrix('Stoffdaten.xlsx', 'Sheet', sheet{prop});
        if mix == 1
            Theta = omegamat(2:17,2:5);                                      % temperature / °C
            omega = omegamat(23:38,2:5);                                       % viscosity of reaction mixtures / mPa s
        elseif mix == 2
            Theta = omegamat(nrm+1,10:13);                                    % temperature / °C
            omega = omegamat(nrm+22, 10:13);                                    % viscosity of artificial mixtures / mPa s
        end
    end
    
    %% fitting paramters
    PropertiesFit

% end