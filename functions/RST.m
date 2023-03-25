function [gamma18, gamma0, gammaH2, gE] = RST(const, Tr, x18, x0)
    % [gamma18, gamma0, gammaH2, gE] = RST(const, Tr, x18, [x0])
    %
    % activity coefficient of H18-DBT (-)
    %                         H0-DBT  (-)
    %                         H2      (-)
    % molar excess Gibbs energy of liquid mixture, containing H2 (kJ mol^-1)
    % 
    % Tr (K)            scalar / row-vector
    % x18 / x0 (-)      scalar / vector / matrix
    %
    % see DOI: 10.1002/aic.690070430
    %          10.1021/acs.jced.5b00789
    %    ISBN: 0-07-0011682-2 (pp. 8.43 f (& 8.32 f))
    %          0-13-977745-8 (pp. 603 ff)
    
    if nargin <= 3
        x0 = 1 - x18;                                                       % molar fraction of H0-DBT in liquid phase (-)
    else
        if ~isequal(size(x18), size(x0))
            error('size(x18) ~= size(x0)')
        end
    end
    xH2 = 1 - x18 - x0;                                                     % molar fraction of H2 in liquid phase (-)
    
    nTheta = numel(Tr);
    if ~isrow(Tr)
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(x18), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(x18) && isrow(x18)
            x18 = x18';
            x0 = x0';
            xH2 = xH2';
        elseif ~isvector(x18) && size(x18,2) ~= nTheta
            error(['If x18 & x0 are matrices, size(x18,2) & size(x0,2) ', ...
                'have to be size(Tr,2)']) 
        end
    end
    
    %% Regular Solution Theory
    % see DOI: 10.1002/aic.690070430
    %          10.1021/acs.jced.5b00789
    %    ISBN: 0-07-0011682-2 (pp. 8.43 f (& 8.32 f))
    %          0-13-977745-8 (pp. 603 ff)
    rho = getRho(const, [], const.T0, [1;0]);                               % density of H18-DBT & H0-DBT (g l^-1)
    v18 = const.M18/rho(1);                                                 % molar volume of H18-DBT @ T0 (l mol^-1)
    v0 = const.M0/rho(2);                                                   % molar volume of H0-DBT @ T0 (l mol^-1)
    vH2 = 31/1000;                                                          % molar volume of H2 in liquid phase @ T0 (l mol^-1)
    
    delta18 = 389.7345;                                                     % solubility paramter of H18-DBT @ T0 (J l^-1)^0.5
    delta0 = 454.5652;                                                      % solubility paramter of H0-DBT @ T0 (J l^-1)^0.5
    deltaH2 = 3.25*sqrt(4.184*1000);                                        % solubility paramter of H2 @ T0 (J l^-1)^0.5
    
    phi18 = x18*v18./(x18*v18 + x0*v0 + xH2*vH2);                           % volume fraction of H18-DBT (-)
    phi0 = x0*v0./(x18*v18 + x0*v0 + xH2*vH2);                              % volume fraction of H0-DBT (-)
    phiH2 = 1 - phi18 - phi0;                                               % volume fraction of H2 (-)
    
    delta = phi18*delta18 + phi0*delta0 + phiH2*deltaH2;                    % solubility paramter of liquid mixture, containing H2 (J l^-1)^0.5
    
    gamma18 = exp(v18*(delta18-delta).^2./(const.R*Tr));                    % activity coefficient of H18-DBT (-)
    gamma0 = exp(v0*(delta0-delta).^2./(const.R*Tr));                       % activity coefficient of H0-DBT (-)
    gammaH2 = exp(vH2*(deltaH2-delta).^2./(const.R*Tr));                    % activity coefficient of H2 (-)

    if nargout == 4
        gE = (const.R*Tr.*(x18.*log(gamma18) + x0.*log(gamma0) + ...        % molar excess Gibbs energy of liquid mixture, containing H2 (kJ mol^-1)
            xH2.*log(gammaH2)))/1000;
    end
end