function [lambdaPe, lambdaPed, lambdaPer] = getLambdaPe(const, p, pr, Tr, HG)
    % [lambdaPe, lambdaPed, lambdaPer] = getLambdaPe(const, p, pr, Tr, HG)
    %
    % effective heat conductivity of porous pellet (catalyst) (W m^-1 K^-1)
    %                                active layer             (W m^-1 K^-1)      
    %                                support material         (W m^-1 K^-1) 
    %
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %       
    % see DOI: 10.3390/ma7117173
    %          10.1557/jmr.2013.179
    %          10.1039/c5cp04745e
    % see ISBN: 3-527-31000-2 (p. 49)
    % see URL: http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/al2o3.html 
    
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
    elseif ~isequal(size(HG), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end

    if ~isfield(p, 'dpore')
        p.dpore = const.dpore;
    end
                                                                     
    lambda = getLambda(const, pr, Tr, HG);                                  % heat conductivity of mixture, wo H2 (W m^-1 K^-1)  
    [rhoAlox, ~, lambdaAlox] = getAlox(const, Tr);                          % solid density (g l^-1) and heat conductivity of alumina (W m^-1 K^-1)
    rhoPt = 21557.19 - 0.5675783*Tr - 1.7525e-05*Tr.^2 - ....               % solid density of Pt (g l^-1)
                3.171806e-08*Tr.^3 + 4.698968e-12*Tr.^4;
    lambdaPt = 73.99627 - 0.01557887*Tr + 2.646931e-05*Tr.^2 - ...          % heat conductivity of Pt (W m^-1 K^-1)
        6.133801e-09*Tr.^3;
    
    if isfield(p,'phiPt')                                                   % volume fraction of catalyst with Pt
        phiPt = p.phiPt;    
    else
        phiPt = const.phiPt;
    end
    if isfield(p,'wPtL')                                                    % mass fraction of Pt regarding to mass of active layer     
        wPtL = p.wPtL;          
    elseif isfield(p,'wPtc')
        wPtL = p.wPtc./phiPt;
    else
        wPtL = const.wPtc./phiPt;
    end
    
    phiPtL = wPtL./rhoPt./(wPtL./rhoPt + (1-wPtL)./rhoAlox);                % volume fraction of Pt regarding to volume of active layer     

    lambdaL = (1-phiPtL).*lambdaAlox + phiPtL.*lambdaPt;                    % heat conductivity of active layer (W m^-1 K^-1)
    lambdaL = repmat(lambdaL,size(HG,1),1);
    
    switch p.lambdaPed                                                      % effective heat conductivity of active layer (catalyst) (W m^-1 K^-1)
        case 1
            %1 - Parallel slabs / rule of mixtures (upper bound)
            lambdaPed = lambdaL.*(1-const.epsPe) + const.epsPe.*lambda; 
        case 2
            %2 - Serial slabs (lower bound)
            lambdaPed = (const.epsPe./lambda + (1-const.epsPe)./lambdaL).^-1;
        case 3
            %3 - Harriot
            delta = 2.*(0.64*3.*(HG.*const.M18 + ...                        % mean free path of liquid = mean molecular spacing (nm)
                (1-HG).*const.M0)/1000./... 
                (4*pi.*getRho(const, pr, Tr, HG).*const.NA)).^(1/3)*1e09; 
            d = 0.1*(p.dpore*1e+09);
            Kn = delta./d;
            Delta = d./(1+d);
            chi = (lambda./(1+2*Kn))./lambdaL;
            lambdaPed = lambdaL.*(Delta.^2 + chi.*(1-Delta).^2 + ...
                (1 + d + 1./(chi.*Delta)).^-1);
    end

    lambdaPer = repmat(lambdaAlox,size(HG,1),1);                            % effective heat conductivity of support material (W m^-1 K^-1)    
    switch p.lambdaPe                                                       % effective heat conductivity of porous pellet (catalyst) (W m^-1 K^-1) 
        case 1
            %1 - Parallel slabs / rule of mixtures (upper bound)
            lambdaPe = lambdaPer.*(1-phiPt) + phiPt.*lambdaPed;             
        case 2
            %2 - Serial slabs (lower bound)
            lambdaPe = ((1-phiPt)./lambdaPer + phiPt./lambdaPed).^-1;
    end
end