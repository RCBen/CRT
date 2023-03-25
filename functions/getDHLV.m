function [dHLV18, dHLV0] = getDHLV(const, Tr)
    % [dHLV18, dHLV0] = getDHLV(const, Tr)
    %
    % heat of vaporization of H18-DBT (kJ mol^-1)
    %                         H0-DBT  (kJ mol^-1)
    %
    % Tr (K)            scalar / vector / matrix
    %
    % see DOI: 10.1002/aic.690190120
 
    %% estimation from the law of corresponding states 
    % see DOI: 10.1002/aic.690190120
    Tr18 = Tr/const.Tc18;                                                   % reduced temperature of H18-DBT (-)
    Hr180 = 7.0837 - 3.1933*Tr18;
    Hr180(Tr18>0.75) = 1.481449 + 36.09175*Tr18(Tr18>0.75) - ...
        101.1801*Tr18(Tr18>0.75).^2 + 113.0982*Tr18(Tr18>0.75).^3 - ...
        46.23935*Tr18(Tr18>0.75).^4 - 0.0090026./(1.01 - Tr18(Tr18>0.75));
    Hr181 = 14.602 - 11.9866*Tr18;
    Hr181(Tr18>0.75) = 24.83238 + 51.88302*Tr18(Tr18>0.75) - ...
        322.8546*Tr18(Tr18>0.75).^2 + 423.1277*Tr18(Tr18>0.75).^3 - ...
        173.8755*Tr18(Tr18>0.75).^4 - 0.0045311./(1.01 - Tr18(Tr18>0.75));
    dHLV18 = const.R*const.Tc18*(Hr180 + const.w18*Hr181)/1000;             % heat of vaporization of H18-DBT (kJ mol^-1)
    
    Tr0 = Tr/const.Tc0;                                                     % reduced temperature of H0-DBT (-)
    Hr00 = 7.0837 - 3.1933*Tr0;
    Hr00(Tr0>0.75) = 1.481449 + 36.09175*Tr0(Tr0>0.75) - ...
        101.1801*Tr0(Tr0>0.75).^2 + 113.0982*Tr0(Tr0>0.75).^3 - ...
        46.23935*Tr0(Tr0>0.75).^4 - 0.0090026./(1.01 - Tr0(Tr0>0.75));
    Hr01 = 14.602 - 11.9866*Tr0;
    Hr01(Tr0>0.75) = 24.83238 + 51.88302*Tr0(Tr0 > 0.75) - ...
        322.8546*Tr0(Tr0>0.75).^2 + 423.1277*Tr0(Tr0>0.75).^3 - ...
        173.8755*Tr0(Tr0>0.75).^4 - 0.0045311./(1.01 - Tr0(Tr0>0.75));
    dHLV0 = const.R*const.Tc0*(Hr00 + const.w0*Hr01)/1000;                  % heat of vaporization of H0-DBT (kJ mol^-1)
    
    %% estimation from the law of corresponding states
    % see ISBN: 0-07-0011682-2 (pp. 7.16 ff)
%     Tr18 = Tr/const.Tc18;                                                  % reduced temperature of H18-DBT (-)
%     tau18 = 1 - Tr18;
%     dHLV182 = const.R*const.Tc18*...                                       % heat of vaporization of H18-DBT (kJ mol^-1)
%         (7.08*tau18.^0.354 + 10.95*const.w18*tau18.^0.456)/1000;
%     
%     Tr0 = Tr/const.Tc0;                                                    % reduced temperature of H0-DBT (-)
%     tau0 = 1 - Tr0;
%     dHLV02 = const.R*const.Tc0*...                                         % heat of vaporization of H0-DBT (kJ mol^-1)
%         (7.08*tau0.^0.354 + 10.95*const.w0*tau0.^0.456)/1000;

    %% estimation from vapor pressure equations (Clausius-Clapeyron)
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 154 ff) 
    %          10.1021/ja0048634 (pp. 7.14 ff) 
%     [pLV18, pLV0] = getPLV(const, Tr);
%     [rho, ~, ~, rhoG] = getRho(const, pr, Tr, [1;0]);
%     v = [const.M18; const.M0]./rho;                                       
%     vG = [const.M18; const.M0]./rhoG;
%    
%     Tr18 = Tr/const.Tc18;                                                  % reduced temperature of H18-DBT (-)
%     tau18 = 1 - Tr18;
%     psi18 = 5.97616 + 1.29874.*tau18.^0.5.*(0.5.*tau18 - 1.5) - ...
%         0.60394.*tau18.^1.5.*(1.5.*tau18 - 2.5) - ...
%         1.06841.*tau18.^4.*(4.*tau18 - 5) + ...
%         const.w18.*(5.03365 + 1.11505.*tau18.^0.5.*(0.5.*tau18 - 1.5) - ...
%         5.41217.*tau18.^1.5.*(1.5.*tau18 - 2.5) - ...
%         7.46628.*tau18.^4.*(4.*tau18 - 5)) + ...
%         const.w18.^2.*(0.64771 + 2.41539.*tau18.^0.5.*(0.5.*tau18 - 1.5) - ...
%         4.26979.*tau18.^1.5.*(1.5.*tau18 - 2.5) + ...
%         3.25259.*tau18.^4.*(4.*tau18 - 5));
%     dHLV18 = (psi18.*(pLV18*1e+05).*((vG(1,:) - v(1,:))/1000)./Tr18)./1000;
%     
%     Tr0 = Tr/const.Tc0;                                                    % reduced temperature of H0-DBT (-)
%     tau0 = 1 - Tr0;
%     psi0 = 5.97616 + 1.29874.*tau0.^0.5.*(0.5.*tau0 - 1.5) - ...
%         0.60394.*tau0.^1.5.*(1.5.*tau0 - 2.5) - ...
%         1.06841.*tau0.^4.*(4.*tau0 - 5) + ...
%         const.w0.*(5.03365 + 1.11505.*tau0.^0.5.*(0.5.*tau0 - 1.5) - ...
%         5.41217.*tau0.^1.5.*(1.5.*tau0 - 2.5) - ...
%         7.46628.*tau0.^4.*(4.*tau0 - 5)) + ...
%         const.w0.^2.*(0.64771 + 2.41539.*tau0.^0.5.*(0.5.*tau0 - 1.5) - ...
%         4.26979.*tau0.^1.5.*(1.5.*tau0 - 2.5) + ...
%         3.25259.*tau0.^4.*(4.*tau0 - 5));
%     dHLV0 = (psi0.*(pLV0*1e05).*((vG(2,:) - v(2,:))/1000)./Tr0)./1000;
end