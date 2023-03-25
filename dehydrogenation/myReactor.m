function p = myReactor(const, p)
    % p = myReactor(const, p)
    % 
    % parameters of reactor PD 2.0
    
    p.Vc = p.mc./const.rhoBed*1000;                                         % volume of catalyst in bed (ml)
%     p.epsBed = const.epsBed;                                              % mean bed porosity (-)
    p.NPe = round((1-const.epsBed)*(p.Vc/1000)./const.V);                   % number of pellets in bed (-)
    
    %% reactor
    LR = 0.2200;                                                            % length of reactor (m)
    BR = 0.1350;                                                            % width of reactor (m)
    nc = 9;                                                                 % number of chambers (-)
    p.BK = 0.0200;                                                          % width of chamber (m)
    p.BRxz = p.BK;                                                          % width of 2Dxz model (m)
    ndam = nc-1;                                                            % number of dams (-)
    Ldam = BR-p.BK;                                                         % length of dam (m)
    Bdam = 0.0050;                                                          % width of dam (m)
    Bplate = 0.0020;                                                        % width of plate (m)
    r1 = 0.0020;
    r2 = 0.0060;
    Axy = LR*BR - 2*p.BK*(p.BK + Bplate) - ndam*Bdam*Ldam + ...             % surface area of xy-plane (m^2)
        2*ndam*(1-1/4*pi)*(r1^2-r2^2);
    if ~isfield(p, 'HR')
        p.HR = p.Vc*10^-6/Axy;                                              % filling height of reactor (m)
    end
    VR = Axy*p.HR;                                                          % reaction volume (contains catalyst) (m^3)
    p.V0 = VR*1e06 - p.Vc*(1-const.epsBed);                                 % initial volume of liquid phase (ml)
    p.epsBed = p.V0/(VR*1e06);                                              % mean bed porosity (-)
    
    p.Ayz = p.BRxz*p.HR;                                                    % cross sectional area in yz-plane (m^2)
    p.DhRxz = 4*p.Ayz/(2*(p.BRxz + p.HR));                                  % hydraulic diameter of 2Dxz model (m)   
    p.LRxz = Axy/p.BRxz;                                                    % length of reactor for 2Dxz model (m)

    ORxz = 2*(p.LRxz*p.HR);                                                 % heat transfer surface area for 2Dxz model (m^2)    
    p.SvRxz = ORxz/VR;                                                      % heat transfer surface area to volume ratio for 2Dxz model (m^-1)

%     OR = (2*((LR-p.BK) + (BR-(p.BK + Bplate))) + ...                        % heat transfer surface area for 3D model (m^2)  
%         2*(ndam*Ldam - (p.BK + Bplate)) - 2*ndam*(2-1/2*pi)*(r1+r2))*p.HR;     
%     p.SvR = OR/VR;                                                          % heat transfer surface area to volume ratio for 3D model (m^-1)  
    
    %% thermocouples
    p.xF = (BR/2 - (p.BRxz + Bplate)) + (0:1:ndam)*((BR + Bdam) + ...       % x-position of thermocouples in F-plane in xz-model (m)
        2*(1-1/4*pi)*(r1^2-r2^2)./p.BRxz);
    p.dxF = [0, p.xF(1) + (1/2 + (0:1:ndam-1))*((BR + Bdam) + ...           % thermocouple zones in F-plane in xz-model (m)
        2*(1-1/4*pi)*(r1^2-r2^2)./p.BRxz)];
    
    xE = sort(reshape(p.xF(1:2:end) + [0; ...
        (BR-p.BRxz)/2 + (1-1/4*pi)*(r1^2-r2^2)./p.BRxz; ...
        -(BR-p.BRxz)/2], 1, 15));
    p.xE = xE(2:end-1);                                                     % x-position of thermocouples in E-plane in xz-model (m)    
    p.zE = [0.0045 0.0135];                                                 % z-position of thermocouples in E-plane in xz-model (m)
end