function p = myCoeff(p)
    % p = myCoeff(p)
    %
    % cestimated coefficients for reaction rate (fitted to experimental 
    % data with myFitLSQnl (tr = 0 - 3600 s, pr = 1 - 2 bar, Tr = 280 - 310 °C)
    % 
    % p.rate (dependency)       0 partial density (rho)
    %                           1 molar concentration (c)
    % p.i(nternal) r(esistance) 0 off
    %                           1 on

    if ~isfield(p,'subdir')
        mydir = pwd;
        idcs = strfind(mydir,'\');
        p.subdir = mydir(1:idcs(end)-1);
    end
    pathname = strcat(p.subdir,'\general\KINETIKfitting\RCB');
    S = load(fullfile(pathname, 'myFitLSQnl.mat'));
    
    name = strcat('ir', num2str(p.ir));
    switch p.rate
        case 0
            conc = 'rho';
        case 1
            conc = 'c';
        case 2
            conc = 'w';    
    end
    
    % check whether p.DPe and lambdaPe are the same in myFitLSQnl as in the argument p
    tPe = isequal(p.DPe, S.(conc).(name).p.DPe);
    tlambdaPe = isequal(p.lambdaPe, S.(conc).(name).p.lambdaPe);
    
    if p.ir == 1 && (tPe == 0 || tlambdaPe == 0)
%         warning('Definition of p.DPe and/or p.lambdaPe are not the same')
        msgfig = msgbox('Definition of p.DPe and/or p.lambdaPe are not the same!','Warning','modal');
        uiwait(msgfig)
        disp('Program execution resumed.');
    end

    p.p_est = S.(conc).(name).mdl.Coefficients.Estimate;
    p.Tm = S.(conc).(name).p.Tm(1);
end