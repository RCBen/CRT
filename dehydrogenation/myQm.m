function [Qm, sigma18, sigma0, sigmaH2] = myQm(const, p, pr, Tr, HG)
    % [Qm, sigma18, sigma0, sigmaH2] = myQm(const, p, pr, Tr, HG)
    %
    % mass source (g l_bed^-1 s^-1)
    % species rate for consumption of H18-DBT (g l_bed^-1 s^-1)
    %                  production     H0-DBT  (g l_bed^-1 s^-1) 
    %                  production     H2      (g l_bed^-1 s^-1) 
    % 
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    
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
        elseif ~isvector(HG) && size(HG,2) ~= nTheta && ~isscalar(Tr)
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end

    nges = 1;                                                               % total molar amount (arbitrary) (mol) 

    n18 = HG.*nges;                                                         % molar amount of H18-DBT in liquid phase (mol)
    n0 = nges-n18;                                                          % molar amount of H0-DBT in liquid phase (mol)
%     if isfield(p, 'LRxz')
%         f = p.epsBed./const.epsBed;                                         % conversion factor for continuous operation (-)
%         f = (1-const.epsBed)./(1-p.epsBed);                                 % conversion factor for continuous operation (-)
%        f = (1-p.epsBed)./(1-0.988)*0.0001/p.xPt;
%          f = p.xPt./(1-mean(p.epsBed))/(0.0001/(1-0.988));
%          f = p.xPt./(1-mean(p.epsBed))/(0.0005/(1-0.939));
%         f = p.xPt./(1-mean(p.epsBed))/0.007908;
%     else
         f = 1;                                                              % conversion factor for batch operation (-)
%     end    
    reff4 = f.*myReactionRate(const, p, pr, Tr, n18, n0);                   % reaction rate for dehydrogenation of H18-DBT (mol l_bed^-1 s^-1)

    sigma18 = const.nu(1).*const.M18.*reff4;                                % species rate for consumption of H18-DBT (g l_bed^-1 s^-1)
    sigma0 = const.nu(2).*const.M0.*reff4;                                  % species rate for production of H0-DBT (g l_bed^-1 s^-1)
    sigmaH2 = const.nu(3).*const.MH2.*reff4;                                % species rate for production of H2 (g l_bed^-1 s^-1)
    
    Qm = sigma18 + sigma0;                                                  % mass source (g l_bed^-1 s^-1)
end