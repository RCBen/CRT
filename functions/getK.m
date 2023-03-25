function K = getK(const, pr, Tr)
    % K = getK(const, pr, Tr)
    %
    % equilibrium constant for dehydrogenation of H18-DBT (-)
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
  
    %% Gibbs-Helmholtz & van't Hoff equation
    syms p T                                                        
    [h18, h0, hH2] = getH(const, p, T);   
    dlnKFun = matlabFunction(int(1000*sum([h18; h0; hH2].*const.nu)./ ...   % logarithmic change of equilibrium constant (-)
        (const.R*T^2), T, const.T0, T), 'Vars', [p T]);
    
    K = exp(const.lnK0 + dlnKFun(pr, Tr));                                  % equilibrium constant for dehydrogenation of H18-DBT (-)
end