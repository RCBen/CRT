function [pLV18, pLV0, fPH2] = getPLV(const, Tr)
    % [pLV18, pLV0, fPH2] = getPLV(const, Tr)
    %
    % vapour pressure of H18-DBT (bar)
    %                    H0-DBT  (bar)
    % Prausnitz-fugacity of H2 (bar)
    %
    % Tr (K)            scalar / vector / matrix
    % 
    % see DOI: 10.1002/aic.690070430
    %    ISBN: 0-07-0011682-2 (p. 7.7 f)
    
    %% Ambrosi-Walton corresponding states method
    % see ISBN: 0-07-0011682-2 (pp. 7.7 f)
%     Tr18 = Tr/const.Tc18;                                                   % reduced temperature of H18-DBT (-)
%     tau18 = 1 - Tr18;
%     f180 = (-5.97616*tau18 + 1.29874*tau18.^1.5 - ...
%         0.60394*tau18.^2.5 - 1.06841*tau18.^5)./Tr18;
%     f181 = (-5.03365*tau18 + 1.11505*tau18.^1.5 - ...
%         5.41217*tau18.^2.5 - 7.46628*tau18.^5)./Tr18;
%     f182 = (-0.64771*tau18 + 2.41539*tau18.^1.5 - ...
%         4.26979*tau18.^2.5 + 3.25259*tau18.^5)./Tr18;
%     pLV18 = const.pc18*exp(f180 + const.w18*f181 + const.w18^2*f182);       % vapour pressure of H18-DBT (bar) 
%     
%     Tr0 = Tr/const.Tc0;                                                     % reduced temperature of H0-DBT (-)
%     tau0 = 1 - Tr0;
%     f00 = (-5.97616*tau0 + 1.29874*tau0.^1.5 - ...
%         0.60394*tau0.^2.5 - 1.06841*tau0.^5)./Tr0;
%     f01 = (-5.03365*tau0 + 1.11505*tau0.^1.5 - ...
%         5.41217*tau0.^2.5 - 7.46628*tau0.^5)./Tr0;
%     f02 = (-0.64771*tau0 + 2.41539*tau0.^1.5 - ...
%         4.26979*tau0.^2.5 + 3.25259*tau0.^5)./Tr0;
%     pLV0 = const.pc0*exp(f00 + const.w0*f01 + const.w0^2*f02);              % vapour pressure of H0-DBT (bar)

    %% experimental data (t <= 250 °C)
%     pLV18 = exp(-67.6068 + 2.3121e-01.*Tr + 4.1100.*1 - ...                 % vapour pressure of H18-DBT (bar) 
%         8.3559e-03.*Tr.*1 - 1.8121e-04.*Tr.^2 + 6.1075e-01.*1.^2)/1000; 
%     pLV18 = exp(238.5708 - 24686.4350./Tr - 29.9036*log(Tr))/1000;
    pLV18 = exp((301.18 - 105895./Tr - 44.338*log(Tr./298.15))/const.R)* ...% TVT
        1e-05; 
%     
%     pLV0 = exp(-67.6068 + 2.3121e-01.*Tr + 4.1100.*0 - ...                  % vapour pressure of H0-DBT (bar) 
%         8.3559e-03.*Tr.*0 - 1.8121e-04.*Tr.^2 + 6.1075e-01.*0.^2)/1000;
%     pLV0 = exp(215.2933 - 24963.6550./Tr - 26.2078*log(Tr))/1000;
    pLV0 = exp((559.01 - 197920./Tr - 201.61*log(Tr./298.15))/const.R)* ... % TVT
        1e-05; 

    %% exponentiational fit to Prausnitz-fugacity data
    % see DOI: 10.1002/aic.690070430
    if nargout == 3                   
        fPH2 = (2.2139e+08*Tr.^-2.2046)*1.01325;                            % Prausnitz-fugacity of H2 (bar)
    end
end