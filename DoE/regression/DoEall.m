    % module            1   PDBoden V1
    %                   2   PDBoden V2
    % result            1   HG      (hydrogenation grade)
    %                   2   QH2     (hydrogen volume flow rate)
    %                   3   PthH2   (thermal power output due to QH2)
    %                   4   Etar    (energetic efficiency of dehydrogenation)
    %                   5   Eta     (energetic efficiency wo pre-heating)
    %                   6   Etaw0   (energetic efficiency with pre-heating)
    %                   7   LB      (change in light boiling products)
    %                   8   HB      (change in high boiling products)
    %                   9   RTD     (residence time distribution)

[mdl.M1.HG, Error.M1.HG] = DoE(1, 1);
[mdl.M1.QH2, Error.M1.QH2] = DoE(1, 2);
[mdl.M1.PtH2, Error.M1.PtH2] = DoE(1, 3);
[mdl.M1.Etar, Error.M1.Etar] = DoE(1, 4);
[mdl.M1.Eta, Error.M1.Eta] = DoE(1, 5);
[mdl.M1.Etaw0, Error.M1.Etaw0] = DoE(1, 6);
[mdl.M1.LB, Error.M1.LB] = DoE(1, 7);
[mdl.M1.HB, Error.M1.HB] = DoE(1, 8);
[mdl.M1.RTD, Error.M1.RTD] = DoE(1, 9);

[mdl.M2.HG, Error.M2.HG] = DoE(2, 1);
[mdl.M2.QH2, Error.M2.QH2] = DoE(2, 2);
[mdl.M2.PtH2, Error.M2.PtH2] = DoE(2, 3);
[mdl.M2.Etar, Error.M2.Etar] = DoE(2, 4);
[mdl.M2.Eta, Error.M2.Eta] = DoE(2, 5);
[mdl.M2.Etaw0, Error.M2.Etaw0] = DoE(2, 6);
%[mdl.M2.LB, Error.M2.LB] = DoE(2, 7);
%[mdl.M2.HB, Error.M2.HB] = DoE(2, 8);
%[mdl.M2.RTD, Error.M2.RTD] = DoE(2, 9);

save('DoE.mat')