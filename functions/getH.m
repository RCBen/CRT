function [h18, h0, hH2] = getH(const, pr, Tr)
    % [h18, h0, hH2] = getH(const, pr, Tr)
    %
    % molar enthalpy of H18-DBT (kJ mol^-1)
    %                   H0-DBT  (kJ mol^-1)
    %                   H2      (kJ mol^-1)
    % 
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    
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

    %% Kirchhoff's law of thermochemistry
    syms p T                                                    
    [cp, cpH2Id] = getCp(const, [], T, [1;0]);                              % isobaric heat capacity of H18-DBT & H0-DBT & ideal gas H2 (kJ mol^-1 K^-1)
    rho = getRho(const, [], T, [1;0]);                                      % density of H18-DBT & H0-DBT (g l^-1)
    beta = -1./rho.*diff(rho,T);                                            % thermal expansion coefficient (K^-1)
%     betaH2 = 1./T;
    
    assume(p <= 10)
    assumeAlso(T >= 0)
    HFun = matlabFunction([const.dHb18; const.dHb0; const.dHbH2] + ...
        int([cp; cpH2Id], const.T0, T) + ...
        int([1./rho.*(1 - beta.*T); 0].* ...
        [const.M18; const.M0; const.MH2], const.p0, p)*1e-06*1e+05, ...
        'Vars', [p T]);
    H = HFun(pr, Tr);
    
    h18 = real(H(1,:));                                                     % molar enthalpy of H18-DBT (kJ mol^-1)
    h0 = real(H(2,:));                                                      % molar enthalpy of H0-DBT (kJ mol^-1)
    hH2 = real(H(3,:));                                                     % molar enthalpy of H2 (kJ mol^-1)
end