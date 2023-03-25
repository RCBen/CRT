function [cp, cpH2Id, cpId] = getCp(const, pr, Tr, HG)
    % [cp, cpH2Id, cpId] = getCp(const, pr, Tr, HG)
    %
    % isobaric heat capacity of liquid mixture, wo H2 (kJ mol^-1 K^-1)
    %                           ideal gas H2          (kJ mol^-1 K^-1)
    %                           ideal gas mixture     (kJ mol^-1 K^-1)
    %
    % pr (bar)          scalar
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 358 & 405)  
    %    ISBN: 0-07-0011682-2 (pp. 3.8 f & C.6 & C.9)
    
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
    
    %% experimental data (t = 30, 70, 110, 150, 180, 210, 230, 250 °C)
    M = HG*const.M18 + (1-HG)*const.M0;                                     % molar mass of liquid mixture, wo H2 (g mol^-1)
    cp = (7.4984e-01 + 2.6880e-03.*Tr - 1.3226e-01.*HG + ...                % isobaric heat capacity of liquid mixture, wo H2 (kJ mol^-1 K^-1)        
        1.1643e-03.*Tr.*HG).*M/1000;    
    
    %% polynomial fit to VDI data
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 358 (eq. 10) & 405)
    cpH2Id = 4.9852e-03 + 2.3150e-04.*Tr - 8.8784e-07.*Tr.^2 + ...          % ideal gas heat capacity of H2 @ constant pressure (kJ mol^-1 K^-1)
        1.7610e-09.*Tr.^3 - 1.9310e-12.*Tr.^4 + 1.12695e-15.*Tr.^5 - ...
        2.7563e-19.*Tr.^6;
    
    %% group contribution method of Constantinou & Gani
    % see ISBN: 0-07-0011682-2 (pp. 3.8 f & C.6 & C.9)   
    if nargout == 3
        theta = (Tr - 298)/700;                                             % reduced temperature (-)
        
        A18 = 1*35.1152 + 15*22.6346 + 5*8.9272 + 3*-3.9682;
        B18 = 1*39.5923 + 15*45.0933 + 5*59.9786 + 3*17.7889;
        C18 = 1*-9.9232 + 15*-15.7033 + 5*-29.5143 + 3*-3.3639;

        cp18Id = ((A18 - 19.7779) + (B18 + 22.5981).*theta + ...            % ideal gas heat capacity of H18-DBT @ constant pressure (kJ mol^-1 K^-1)             
            (C18 - 10.7983).*theta.^2)/1000;   

        A0 =  13*16.3794 + 2*10.4283 + 1*42.8569 + 2*32.8206;
        B0 =  13*32.7433 + 2*25.3634 + 1*65.6464 + 2*70.4153;
        C0 = 13*-13.1692 + 2*-12.7283 + 1*-21.0670 + 2*-28.9361;

        cp0Id = ((A0 - 19.7779) + (B0 + 22.5981).*theta + ...               % ideal gas heat capacity of H0-DBT @ constant pressure (kJ mol^-1 K^-1)
            (C0 - 10.7983).*theta.^2)/1000;

        %% ideal gas law
        if size(HG, 1) == 2 && all(HG == [1;0])
            cpId = [cp18Id; cp0Id];                                         % ideal gas heat capacity of H18-DBT & H0-DBT @ constant pressure (kJ mol^-1 K^-1) 
        else
            vle = VLE(const, pr, Tr, HG);           
            y18 = vle.y18;                                                  % molar fraction of H18-DBT in gas (vapour) phase (-) 
            y0 = vle.y0;                                                    % molar fraction of H0-DBT in gas (vapour) phase (-)
            yH2 = vle.yH2;                                                  % molar fraction of H2 in gas (vapour) phase (-)

            cpId = y18.*cp18Id + y0.*cp0Id + yH2.*cpH2Id;                   % isobaric heat capacity of ideal gas mixture (kJ mol^-1 K^-1)
        end
    end
end