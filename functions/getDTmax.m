function dTmax = getDTmax(const, p, pr, Tr, HG)
    % dTmax = getDTmax(const, p, pr, Tr, HG)
    %
    % maximum temperature difference inside the active layer (K)
    %
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %       
    % see ISBN: 3-527-31000-2 (p. 124)
    
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
    if ~isrow(Tr)
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr))
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end

    dHr = getDHr(const, pr, Tr);                                            % reaction enthalpy for dehydrogenation of H18-DBT (kJ mol^-1)
    DPe = getDPe(const, p, Tr, HG);                                         % binary diffusion coefficient of H18-DBT & H0-DBT in liquid mixture, wo H2 (m^2 s^-1)
    rho = getRho(const, pr, Tr, HG);                                        % density of liquid mixture, wo H2 (g l^-1)
%     lambdaPe = getLambdaPe(const, p, pr, Tr, HG);                           % effective heat conductivity of porous pellet (catalyst) (W m^-1 K^-1)
    [~, lambdaPed] = getLambdaPe(const, p, pr, Tr, HG);                     % effective heat conductivity of active layer (catalyst) (W m^-1 K^-1)
    
    dTmax = 1000*dHr.*DPe.* ...
        (1000*rho./(HG.*const.M18 + (1-HG).*const.M0))./lambdaPed;
end