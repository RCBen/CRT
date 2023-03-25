function dHr = getDHr(const, pr, Tr)
    % dHr = getDHr(const, pr, Tr)
    %
    % reaction enthalpy for dehydrogenation of H18-DBT (kJ mol^-1)
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

    [h18, h0, hH2] = getH(const, pr, Tr);                                   % molar enthalpy of H18-DBT, H0-DBT & H2 (kJ mol^-1)
    dHr = sum([h18; h0; hH2].*const.nu);                                    % reaction enthalpy for dehydrogenation of H18-DBT (kJ mol^-1)
end