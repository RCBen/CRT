function [ZG, phi18V, phi0V, phiH2G, phi18LV, phi0LV, phiH2LV] = ...
    SRK(const, pr, Tr, y18, y0)
    % [ZG, phi18V, phi0V, phiH2G, phi18LV, phi0LV, phiH2LV] = ...
    %     SRK(const, pr, Tr, y18, y0)
    %
    % compressibility factor of real gas mixture (-)
    % fugacity coefficient of H18-DBT (-) @ pr (V/G)
    %                         H0-DBT  (-)
    %                         H2      (-) 
    % fugacity coefficient of H18-DBT (-) @ saturation (LV)
    %                         H0-DBT  (-)
    %                         H2      (-)
    %                         
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % y18 / y0 (-)      scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 150 f & 591 f)
    %          10.1021/i260068a009
    %          10.1021/i260070a022
    %    ISBN: 0-07-0011682-2 (pp. 4.17 ff)
    
    nP = numel(pr);
    nTheta = numel(Tr);
    if ~isrow(pr)
        error('pr has to be a scalar or a row-vector')
    elseif isrow(pr) && (nP ~= nTheta || ~isscalar(pr) && pr(1) ~= pr(2))
        pr = repelem(pr,1,nTheta);
        Tr = repmat(Tr,1,nP);
    end
    if ~isrow(Tr)
        error('Tr has to be a row-vector')
    elseif ~isequal(size(y18), size(y0))
        error('size(y18) ~= size(y0)')
    elseif ~isequal(size(y18), size(Tr))
        if ~isscalar(y18) && isrow(y18)
            y18 = y18';
            y0 = y0';
        elseif ~isvector(y18) && size(y18,2) ~= nTheta
            error(['If y18 & y0 are matrices, size(y18,2) & size(y0,2) ', ...
                'have to be size(Tr,2)']) 
        end
    end
    
    optionsFzero = optimset('Display','off');
    
    [pLV18, pLV0] = getPLV(const, Tr);                                      % vapour pressure of H18-DBT & H0-DBT (bar)
    yH2 = 1 - y18 - y0;                                                     % molar fraction of H2 in gas (vapour) phase (-)
    
    %% Soave-Redlich-Kwong EOS
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 150 f & 591 f)
    %          10.1021/i260068a009
    %          10.1021/i260070a022
    %    ISBN: 0-07-0011682-2 (pp. 4.17 ff)
    a18 = 1/(9*(2^(1/3) - 1))*const.R^2*const.Tc18^2/(const.pc18*1e05)* ...        
        (1 + (0.48508 + 1.55171*const.w18 - 0.15613*const.w18^2)* ...
        (1 - sqrt(Tr/const.Tc18))).^2;
    b18 = (2^(1/3) - 1)/3*const.R*const.Tc18/(const.pc18*1e05);
    A18 = a18.*(pLV18*1e05)./(const.R*Tr).^2;
    B18 = b18*(pLV18*1e05)./(const.R*Tr);  
    
    a0 = 1/(9*(2^(1/3) - 1))*const.R^2*const.Tc0^2/(const.pc0*1e05)* ...
        (1 + (0.48508 + 1.55171*const.w0 - 0.15613*const.w0^2)* ...
        (1 - sqrt(Tr/const.Tc0))).^2;
    b0 = (2^(1/3) - 1)/3*const.R*const.Tc0/(const.pc0*1e05);
    A0 = a0.*(pLV0*1e05)./(const.R*Tr).^2;
    B0 = b0*(pLV0*1e05)./(const.R*Tr); 

    aH2 = 1/(9*(2^(1/3) - 1))*const.R^2*const.TcH2^2/(const.pcH2*1e05)* ...
        (1.202*exp(-0.30228*(Tr/const.TcH2)));
    bH2 = (2^(1/3) - 1)/3*const.R*const.TcH2/(const.pcH2*1e05);
   
    k018 = 0;
    kH218 = 0;
    k180 = 0;
    kH20 = 0;
    k18H2 = 0;
    k0H2 = 0;
    a = y18.*y18.*sqrt(a18.*a18) + y0.*y18.*sqrt(a0.*a18)*(1-k018) + ...
        yH2.*y18.*sqrt(aH2.*a18)*(1-kH218) + ...
        y18.*y0.*sqrt(a18.*a0)*(1-k180) + y0.*y0.*sqrt(a0.*a0) + ...
        yH2.*y0.*sqrt(aH2.*a0)*(1-kH20) + ...
        y18.*yH2.*sqrt(a18.*aH2)*(1-k18H2) + ...
        y0.*yH2.*sqrt(a0.*aH2)*(1-k0H2) + yH2.*yH2.*sqrt(aH2.*aH2);
    b = y18*b18 + y0*b0 + yH2*bH2;
    A = a.*(pr*1e05)./(const.R*Tr).^2;
    B = b.*(pr*1e05)./(const.R*Tr);
    
    Z18 = zeros(size(B18));                                                 % compressibility factor of H18-DBT (-)
    for j = 1:size(B18,2)
        Z18(1,j) = fzero(@(Z18) ...
            1*Z18^3 - ...
            1*Z18^2 + ...
            (A18(1,j) - B18(1,j)^2 - B18(1,j))*Z18 - ...
            (A18(1,j)*B18(1,j)), 1, optionsFzero);
    end
    v18LV = Z18*const.R.*Tr./(pLV18*1e05);                                  % molar volume of H18-DBT in gas (vapour) phase @ pLV18 (saturation) (m^3 mol^-1)
    phi18LV = exp((Z18-1) - log(Z18.*(1-b18./v18LV)) - ...                  % fugacity coefficient of H18-DBT @ pLV18 (saturation) (-)
        a18./(b18*const.R*Tr).*log(1+b18./v18LV)); 
    
    Z0 = zeros(size(B0));                                                   % compressibility factor of H0-DBT (-)
    for j = 1:size(B0,2)
        Z0(1,j) = fzero(@(Z0) ...
            1*Z0^3 - ...
            1*Z0^2 + ...
            (A0(1,j) - B0(1,j)^2 - B0(1,j))*Z0 - ...
            (A0(1,j)*B0(1,j)), 1, optionsFzero);
    end
    v0LV = Z0*const.R.*Tr./(pLV0*1e05);                                     % molar volume of H0-DBT in gas (vapour) phase @ pLV0 (saturation) (m^3 mol^-1)    
    phi0LV = exp((Z0-1) - log(Z0.*(1-b0./v0LV)) - ...                       % fugacity coefficient of H0-DBT @ pLV0 (saturation) (-)
        a0./(b0*const.R*Tr).*log(1+b0./v0LV));
    
    phiH2LV = ones(nTheta);                                                 % fugacity coefficient of H2 @ fPH2 (saturation) (-)
                                                                            % = 1 because already included in Prausnitz-fugacity (fPH2)
                                                                            
    ZG = zeros(size(B));                                                    % compressibility factor of real gas mixture (-)
    for i =1:size(B,1)
        for j = 1:size(B,2)
            ZG(i,j) = fzero(@(Z) ...
                1*Z^3 - ...
                1*Z^2 + ...
                (A(i,j) - B(i,j)^2 - B(i,j))*Z - ...
                (A(i,j)*B(i,j)), 1, optionsFzero);
        end
    end 
    vG = ZG*const.R.*Tr./(pr*1e05);                                         % molar volume of real gas mixture (m^3 mol^-1)
    phi18V = exp(b18./b.*(ZG-1) - log(ZG.*(1-b./vG)) + ...                  % fugacity coefficient of H18-DBT @ pr (real gas mixture) (-)
        1./(b*const.R.*Tr).*(a*b18./b - 2*sqrt(a.*a18)).*log(1+b./vG)); 
    phi0V = exp(b0./b.*(ZG-1) - log(ZG.*(1-b./vG)) + ...                    % fugacity coefficient of H0-DBT @ pr (real gas mixture) (-)
        1./(b*const.R.*Tr).*(a*b0./b - 2*sqrt(a.*a0)).*log(1+b./vG));
    phiH2G = exp(bH2./b.*(ZG-1) - log(ZG.*(1-b./vG)) + ...                  % fugacity coefficient of H2 @ pr (real gas mixture) (-)
        1./(b*const.R.*Tr).*(a*bH2./b - 2*sqrt(a.*aH2)).*log(1+b./vG));
end