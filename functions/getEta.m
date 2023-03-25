function [eta, etaH2Id, etaId, etaG]  = getEta(const, pr, Tr, HG)
    % [eta, etaH2Id, etaId, etaG]  = getEta(const, pr, Tr, HG)
    %
    % dynamic viscosity of liquid mixture, wo H2 (mPas)
    %                      ideal gas H2          (uPas)
    %                      ideal gas mixture     (uPas)
    %                      real gas mixture      (uPas)
    %
    % pr (bar)          scalar
    % Tr (K)            scalar / row-vector
    % HG (-)            scalar / vector / matrix
    %
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 162 & 164 & 165 & 357 & 425)
    %    ISBN: 0-07-0011682-2 (pp. 9.9 & 9.23 & 9.35 f & 9.47)
    
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
    if ~isrow(Tr) && ~isequal(size(HG), size(Tr))
        error('Tr has to be a scalar or a row-vector')
    elseif ~isequal(size(HG), size(Tr)) && ~isscalar(Tr)
        if ~isscalar(HG) && isrow(HG)
            HG = HG';
        elseif ~isvector(HG) && size(HG,2) ~= nTheta
            error('If HG is a matrix, size(HG,2) has to be size(Tr,2)')
        end
    end
      
    %% experimental data (t = 40, 50, 60, 70, 80, 100, 120, 140 °C)
    % see DOI: 10.1007/978-3-642-19981-3 (p. 162)
    % PPDS
%     A = -44.2912 - 617.0057*HG + 733.7414*HG.^2 - 201.0445*HG.^3 + ...
%         703.1453*mod(HG, 2/3) - 1156.9498*mod(HG, 2/3).^2;
%     B = -42.3273 - 617.5978*HG + 736.0528*HG.^2 - 203.8625*HG.^3 + ...
%         711.3167*mod(HG, 2/3) - 1170.2169*mod(HG, 2/3).^2;
%     C = 208.5763 + 274.8677*HG - 748.4242*HG.^2 + 559.3029*HG.^3 - ...
%         161.5070*mod(HG, 2/3) + 251.3101*mod(HG, 2/3).^2;
%     D = 191.7603 + 265.0586*HG - 614.0913*HG.^2 + 425.3200*HG.^3 - ...
%         168.7925*mod(HG, 2/3) + 262.7307*mod(HG, 2/3).^2;
%     E = 9.9991e-03 + 1.1751e-01*HG - 2.9935e-01*HG.^2 + ...
%         2.0789e-01*HG.^3 - 1.3338e-01*mod(HG, 2/3) + ...
%         2.0709e-01*mod(HG, 2/3).^2;
    A = -45.6204 - 596.8185*HG + 703.3793*HG.^2 - 185.9968*HG.^3 + ...
        687.4892*mod(HG, 2/3) - 1142.6758*mod(HG, 2/3).^2;
    B = -43.6920 - 597.5613*HG + 706.2347*HG.^2 - 189.2150*HG.^3 + ...
        695.8433*mod(HG, 2/3) - 1156.1697*mod(HG, 2/3).^2;
    C = 207.0204 + 270.3903*HG - 728.5729*HG.^2 + 544.2787*HG.^3 - ...
        156.2035*mod(HG, 2/3) + 243.6483*mod(HG, 2/3).^2;
    D = 190.4931 + 260.0106*HG - 594.6198*HG.^2 + 410.9372*HG.^3 - ...
        163.2974*mod(HG, 2/3) + 254.9728*mod(HG, 2/3).^2;
    E = 9.9991e-03 + 1.1672e-01*HG - 2.9715e-01*HG.^2 + ...
        2.0640e-01*HG.^3 - 1.3282e-01*mod(HG, 2/3) + ...
        2.0597e-01*mod(HG, 2/3).^2;
    eta = E.*exp(A.*-((Tr - C)./(Tr - D)).^(1/3) + ...                      % dynamic viscosity of liquid mixture, wo H2 (mPas)
        B.*-((Tr - C)./(Tr - D)).^(1/3).*((C - Tr)./(Tr - D)));

    % polynom33
%     eta = exp(36.7814 - ...                                                 % dynamic viscosity of liquid mixture, wo H2 (mPas)
%         2.0623e-01.*Tr + 9.7511.*HG - 3.6930e-02.*Tr.*HG + ...      
%         3.9183e-04.*Tr.^2 + 1.0299.*HG.^2 + ...
%         3.1397e-05.*Tr.^2.*HG + 3.1636e-03.*Tr.*HG.^2 - ...
%         2.5581e-07.*Tr.^3 - 1.8381.*HG.^3);                                                            

    %% VDI data
    % see DOI: 10.1007/978-3-642-19981-3 (pp. 357 (eq. 3) & 425)
    etaH2Id = (0.18024e-05 + 0.27174e-07*Tr - 0.13395e-10*Tr.^2+...         % ideal gas dynamic viscosity of H2 (uPas)
        0.00585e-12*Tr.^3 - 0.00104e-15*Tr.^4)*1e+06;
    
    %% Lucas corresponding states method
    % see DOI: 10.1007/978-3-642-19981-3 (p. 164)
    %    ISBN: 0-07-0011682-2 (p 9.9)
    if nargout >= 3
        Tr18 = Tr/const.Tc18;                                               % reduced temperature of H18-DBT (-)
        mur18 = 52.46*(const.mu18^2*const.pc18*const.Tc18^-2);              % reduced dipole moment of H18-DBT (-)
        Xi18 = 0.176*(const.Tc18*const.M18^-3*const.pc18^-4)^(1/6);         % reduced inverse viscosity of H18-DBT (0.1 uPas)^-1              
        Fp18Id = 1;                                                         % correction term for polarity of H18-DBT
        Fq18Id = 1;                                                         % correction term for quantum effects of H18-DBT    
        eta18Id = ((0.807*Tr18.^0.618 - 0.357*exp(-0.449*Tr18) + ...        % ideal gas dynamic viscosity of H18-DBT (0.1 uPas)
            0.340*exp(-4.058*Tr18) + 0.018)*Fp18Id*Fq18Id/Xi18);

        Tr0 = Tr/const.Tc0;                                                 % reduced temperature of H0-DBT (-)
        mur0 = 52.46*(const.mu0^2*const.pc0*const.Tc0^-2);                  % reduced dipole moment of H0-DBT (-)
        Xi0 = 0.176*(const.Tc0*const.M0^-3*const.pc0^-4)^(1/6);             % reduced inverse viscosity of H0-DBT (0.1 uPas)^-1 
        Fp0Id = 1;                                                          % correction term for polarity of H0-DBT (-)
        Fq0Id = 1;                                                          % correction term for quantum effects of H0-DBT (-)
        eta0Id = ((0.807*Tr0.^0.618 - 0.357*exp(-0.449*Tr0) + ...           % ideal gas dynamic viscosity of H0-DBT (0.1 uPas)
            0.340*exp(-4.058*Tr0) + 0.018)*Fp0Id*Fq0Id/Xi0);

        TrH2 = Tr/const.TcH2q;                                              % reduced temperature of H2 (-)
        murH2 = 52.46*(const.muH2^2*const.pcH2q*const.TcH2q^-2);            % reduced dipole moment of H2 (-)
        FpH2Id = 1;                                                         % correction term for polarity of H2 (-)
        FqH2Id = 1.22*0.76^0.15* ...                                        % correction term for quantum effects of H2 (-)
            (1 + 0.00385*((TrH2 - 12).^2).^(1/const.MH2).*sign(TrH2 - 12));          

        %% method of Reichenberg
        % see ISBN: 0-07-0011682-2 (pp. 9.15 f)
        if nargout == 3 && size(HG, 1) == 2 && all(HG == [1;0])
            etaId = [eta18Id; eta0Id]*0.1;                                  % ideal gas dynamic viscosity of H18-DBT & H0-DBT (uPas)
        else
            vle = VLE(const, pr, Tr, HG);           
            y18 = vle.y18;                                                  % molar fraction of H18-DBT in gas (vapour) phase (-) 
            y0 = vle.y0;                                                    % molar fraction of H0-DBT in gas (vapour) phase (-)
            yH2 = vle.yH2;                                                  % molar fraction of H2 in gas (vapour) phase (-)

            Tr180 = Tr/sqrt(const.Tc18*const.Tc0);
            Tr18H2 = Tr/sqrt(const.Tc18*const.TcH2q);
            Tr0H2 = Tr/sqrt(const.Tc0*const.TcH2q);

            mur180 = sqrt(mur18*mur0);
            mur18H2 = sqrt(mur18*murH2);
            mur0H2 = sqrt(mur0*murH2);

            Fr18 = (Tr18.^3.5 + (10*mur18)^7)./ ...
                (Tr18.^3.5*(1 + (10*mur18)^7));                     
            Fr0 = (Tr0.^3.5 + (10*mur0)^7)./ ... 
                (Tr0.^3.5*(1 + (10*mur0)^7));
            FrH2 = (TrH2.^3.5 + (10*murH2)^7)./ ...
                (TrH2.^3.5*(1 + (10*murH2)^7));

            Fr180 = (Tr180.^3.5 + (10*mur180)^7)./ ...                           
                (Tr180.^3.5*(1 + (10*mur180)^7));
            Fr18H2 = (Tr18H2.^3.5 + (10*mur18H2)^7)./ ...                           
                (Tr18H2.^3.5*(1 + (10*mur18H2)^7));
            Fr0H2= (Tr0H2.^3.5 + (10*mur0H2)^7)./ ...                           
                (Tr0H2.^3.5*(1 + (10*mur0H2)^7));

            U18 = (1 + 0.36*Tr18.*(Tr18 - 1)).^(1/6).*Fr18./sqrt(Tr18);
            U0 = (1 + 0.36*Tr0.*(Tr0 - 1)).^(1/6).*Fr0./sqrt(Tr0);
            UH2 = (1 + 0.36*TrH2.*(TrH2 - 1)).^(1/6).*FrH2./sqrt(TrH2);

            C18 = const.M18^(1/4)./sqrt(eta18Id.*U18);
            C0 = const.M0^(1/4)./sqrt(eta0Id.*U0);
            CH2 = const.MH2^(1/4)./sqrt((etaH2Id*10).*UH2);

            H180 = sqrt(const.M18*const.M0/ ...
                (32*(const.M18 + const.M0)^3))*(C18 + C0).^2.* ...
                ((1 + 0.36*Tr180.*(Tr180-1)).^(1/6).*Fr180./sqrt(Tr180));
            H018 = H180;
            H18H2 = sqrt(const.M18*const.MH2/ ...
                (32*(const.M18 + const.MH2)^3))*(C18 + CH2).^2.* ...
                ((1 + 0.36*Tr18H2.*(Tr18H2-1)).^(1/6).*Fr18H2./sqrt(Tr18H2));
            HH218 = H18H2;
            H0H2 = sqrt(const.M0*const.MH2/ ...
                (32*(const.M0 + const.MH2)^3))*(C0 + CH2).^2.* ...
                ((1 + 0.36*Tr0H2.*(Tr0H2-1)).^(1/6).*Fr0H2./sqrt(Tr0H2));
            HH20 = H0H2;

            K18 = y18.*eta18Id./(y18 + eta18Id.* ...
                (y0.*H180*(3 + (2*const.M0/const.M18)) + ...
                yH2.*H18H2*(3 + (2*const.MH2/const.M18))));
            K0 = y0.*eta0Id./(y0 + eta0Id.* ...
                (y18.*H018*(3 + (2*const.M18/const.M0)) + ...
                yH2.*H0H2*(3 + (2*const.MH2/const.M0))));
            KH2 = yH2.*(etaH2Id*10)./(yH2 + (etaH2Id*10).* ...
                (y18.*HH218*(3 + (2*const.M18/const.MH2)) + ...
                y0.*HH20*(3 + (2*const.M0/const.MH2))));            

            etaId = (K18.*(1 + 2*0 + ((H180.*K0).^2 + ...                   % dynamic viscosity of ideal gas mixture (uPas)
                2*H180.*H18H2.*K0.*KH2 + (H18H2.*KH2).^2)) + ...
                K0.*(1 + 2*(H018.*K18) + ((H018.*K18).^2 + ...
                2*H018.*H0H2.*K18.*KH2 + (H0H2.*KH2).^2)) + ...
                KH2.*(1 + 2*(HH218.*K18 + HH20.*K0) + ((HH218.*K18).^2 + ...
                2*HH218.*HH20.*K18.*K0 + (HH20.*K0).^2)))*0.1;
        end
        
        %% Lucas corresponding states method
        % see DOI: 10.1007/978-3-642-19981-3 (p. 165)
        %    ISBN: 0-07-0011682-2 (pp. 9.47 (& 9.9 & 9.23 & 9.35 f))
        if nargout == 4
            Tc = y18*const.Tc18 + y0*const.Tc0 + yH2*const.TcH2q;           % critical temperature of gas mixture (K)
            Zc = y18*const.Zc18 + y0*const.Zc0 + yH2*const.ZcH2q;           % critical compressibility factor of gas mixture (-)
            vc = y18*const.vc18 + y0*const.vc0 + yH2*const.vcH2q;           % critical molar volume of gas mixture (l mol^-1)
            pc = (const.R.*Tc.*Zc./(vc/1000))*1e-05;                        % critical pressure of gas mixture (bar)

            FpId = y18.*Fp18Id + y0.*Fp0Id + yH2.*FpH2Id;                   % correction term for polarity of ideal gas mixture (-)
            Aq = ones(size(y18));
            Aq(yH2 > 0 & y18 > 0.05 & y18 < 0.7) = ...
                1 - 0.01*(const.M18/const.MH2)^0.87;
            FqId = (y18.*Fq18Id + y0.*Fq0Id + yH2.*FqH2Id).*Aq;             % correction term for quantum effects of ideal gas mixture (-)
            
            a = 1.245e-03./(Tr./Tc).*exp(5.1726*(Tr./Tc).^-0.3286);
            b = a.*(1.6553*(Tr./Tc) - 1.2723);
            c = 0.4489./(Tr./Tc).*exp(3.0578*(Tr./Tc).^-37.7332);
            d = 1.7368./(Tr./Tc).*exp(2.2310*(Tr./Tc).^-7.6351);
            e = 1.3088;
            f = 0.9425*exp(-0.1853*(Tr./Tc).^0.4489);
            Y = 1 + (a.*(pr./pc).^e./ ...
                (b.*(pr./pc).^f + (1 + c.*(pr./pc).^d).^-1));
            FpG = (1 + (FpId - 1).*Y.^-3)./FpId;                            % correction term for polarity of real gas mixture (-)
            FqG = (1 + (FqId - 1).*(Y.^-1 - 0.007*log(Y).^4))./FqId;        % correction term for quantum effects of real gas mixture (-)    
            
            etaG = etaId.*Y.*FpG.*FqG;                                      % dynamic viscosity of real gas mixture (uPas)            
        end
    end
end