function dTmaxAnderson = getDTmaxAnderson(const, p, Tr)
    % dTmaxAnderson = getDTmaxAnderson(const, p, Tr)
    %
    % maximum temperature difference to avoid external and internal heat transfer limitations in pellet (K)
    %
    % Tr (K)            scalar / row-vector
    % 
    % Anderson criteria: 1-k/kB < 5%
    
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

    dTmaxAnderson = Tr - (1./(Tr) - log(1-0.05).*const.R./S.(conc).(name).r.EA).^-1;
end