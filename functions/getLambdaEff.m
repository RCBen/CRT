function lambdaEff = getLambdaEff(const, p, pr, Tr, HG)
    % lambdaEff = getLambdaEff(const, p, pr, Tr, HG)
    %
    % effective heat conductivity of catalyst bed (W m^-1 K^-1)
    %
    % pr (bar)          scalar / row-vector
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %       
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 651 ff)
    
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
                                                                     
    lambda = getLambda(const, pr, Tr, HG);                                  % heat conductivity of mixture, wo H2 (W m^-1 K^-1)
    lambdaPe = getLambdaPe(const, p, pr, Tr, HG);                           % effective heat conductivity of porous pellet (catalyst) (W m^-1 K^-1)
    
    switch p.lambdaEff
        case 1
            %1 - Parallel slabs / rule of mixtures (upper bound)
            lambdaEff = lambdaPe.*(1-p.epsBed) + p.epsBed.*lambda;
        case 2
            %2 - Serial slabs (lower bound)
            lambdaEff = (p.epsBed./lambda + (1-p.epsBed)./lambdaPe).^-1;
        case 3
            %3 - Power law
            lambdaEff = lambda.^p.epsBed.*lambdaPe.^(1-p.epsBed);
        otherwise
            kp = lambdaPe./lambda;
            Cf = 1.25;
            B = Cf.*((1-p.epsBed)./p.epsBed).^(10/9).*1;
            kG = 1;           
            switch p.lambdaEff
                case 4
                    %4 - Zehner/Bauer/Schlünder (primary)
                    krad = 0;
                    phi = 0;
                case 5
                    %5 - Zehner/Bauer/Schlünder (secondary)
                    epsilon = 0.9;
                    krad = 4.*const.sigma./(2/epsilon - 1).*Tr.^3.* ...
                        (2.*const.rd)./lambda;
                    phi = 0.0077;
            end
            N = 1./kG.*(1 + (krad - B.*kG)./kp) - ...
                B.*(1./kG - 1).*(1 + krad./kp);
            kc = 2./N.*(B.*(kp + krad - 1)./(N.^2.*kG.*kp).*...
                log((kp + krad)./(B.*(kG + (1-kG).*(kp+krad)))) + ...
                (B+1)./(2.*B).*(krad./kG - B.*(1 + (1-kG)./kG.*krad)) - ...
                (B-1)./(N.*kG));
            lambdaEff = lambda.*((1 - sqrt(1-p.epsBed)).*p.epsBed.*...
                ((p.epsBed - 1 + kG.^-1).^-1 + krad) + ...
                sqrt(1-p.epsBed).*(phi.*kp + (1-phi).*kc));
    end
end