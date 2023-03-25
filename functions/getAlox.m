function [rhoAlox, cpAlox, lambdaAlox] = getAlox(const, Tr)
    % [rhoAlox, cpAlox, lambdaAlox] = getAlox(const, Tr)
    %
    % solid density          of alumina (g l^-1)
    % isobaric heat capacity of alumina (J kg^-1 K^-1)
    % heat conductivity      of alumina (W m^-1 K^-1)
    %
    % Tr (K)            scalar / vector / matrix
        
%     syms T                                                        
%     alpha = (-0.23036 + 7.0045e-04*T + 5.681*10^-08*T^2)*1e-06;
%     
%     rhoAloxFun = matlabFunction(const.rhoAlox* ...
%         exp(-int(3*alpha, const.T0, T)));
%     rhoAlox = rhoAloxFun(Tr);                                               % solid density of alumina (g l^-1)    
    rhoAlox = const.rhoAlox;                                                % solid density of alumina (g l^-1)
    
    cpAlox = -40.92 + 4.024.*Tr - 5.0048e-03.*Tr.^2 + ...                   % isobaric heat capacity of alumina (J kg^-1 K^-1)
        2.8852e-06.*Tr.^3 - 6.2488e-10.*Tr.^4;
    
    lambdaAlox = 85.868 - 0.22972.*Tr + 2.607e-04.*Tr.^2 - ...              % heat conductivity of alumina (W m^-1 K^-1)
        1.3607e-07.*Tr.^3 + 2.7092e-11.*Tr.^4;
end