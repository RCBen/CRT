function p = myDefinitions(p)
    % p = myDefinitions(p)
    %
    % standardized defintion of below parameters
    
    p.rate = 2;                                                             % rate (dependency)             0 partial density (rho)
                                                                            %                               1 molar concentration (c)
                                                                            %                               2 weight fraction (w)
    p.ir = 0;                                                               % i(nternal) r(esistance)       0 off
                                                                            %                               1 on 
    p.irT = 0;                                                              % ir regarding temperature      0 off
                                                                            %                               1 on       
    p.film = 0;                                                             % film transport limitation     0 off
                                                                            %                               1 on                                                                      
    p.iso = 0;                                                              % iso(thermal)                  0 off
                                                                            %                               1 on
    p.DPe = 3;                                                              % correction factor             1 ratio 1/10
                                                                            %                               2 Dumanski
                                                                            %                               3 Hugo
                                                                            %                               4 Bruggemann
                                                                            %                               5 Maxwell
                                                                            %                               6 Millington/Quirk        
    p.FPe = 0;                                                              % restrictive diffusion         0 off
                                                                            %                               1 on
    p.DEff = 4;                                                             % correction factor (DEff)      1 ratio 1/10
                                                                            %                               2 Dumanski
                                                                            %                               3 Hugo
                                                                            %                               4 Bruggemann
                                                                            %                               5 Maxwell
                                                                            %                               6 Millington/Quirk 
    p.lambdaPed = 2;                                                        % correction factor (lambdaPed) 1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)
                                                                            %                               3 Harriot
    p.lambdaPe = 1;                                                         % correction factor (lambdaPe)  1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)                                                                        
    p.lambdaEff = 1;                                                        % correction factor (lambdaEff) 1 parallel slabs (upper bound)
                                                                            %                               2 serial slabs (lower bound)
                                                                            %                               3 power law
                                                                            %                               4 Zehner/Bauer/Schlünder (primary)
                                                                            %                               5 Zehner/Bauer/Schlünder (secondary)
end
