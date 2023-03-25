function [Q, T, ypred] = Predict(module, result, p) 
    % [Q, T, ypred] = Predict(module, result, p) 
    % 
    % module            1   PDBoden V1
    %                   2   PDBoden V2
    % result            1   HG      (hydrogenation grade)
    %                   2   QH2     (hydrogen volume flow rate)
    %                   3   PthH2   (thermal power output due to QH2)
    %                   4   eta     (energetic efficiency wo pre-heating)
    %                   5   etaw0   (energetic efficiency with pre-heating)
    %                   6   LB      (change in light boiling products)
    %                   7   HB      (change in high boiling products)
    %                   8   RTD     (residence time distribution)
    % p(ressure)        bara
    
    load('DoE')
    if nargin <= 2
        p = 1;
        if nargin <= 1
            result = 3;
            if nargin == 0
                module = 1;
            end
        end
    end    
    plotting = 0;
    
    N = 50;
    Ncontour = 10;
    Qmin = 10;
    Qmax = 20;
    pmin = 1; 
    pmax = 2; 
    Tmin = 280;
    Tmax = 310;
    
    nmin = -1;
    nmax = +1;
    
    mp = (nmax-nmin)/(pmax-pmin);
    tp = -mp*mean([pmax, pmin]);
    if p < pmin || p > pmax
        error('p has to be between pmin and pmax')
    end
    
    Variables = {'QLOHC' 'p' 'T' 'HG' 'QH2' 'PtH2' 'eta' 'etaw0' ...
        'LB' 'HB' 'RTD'};
    Units = {'[ml min^-1]' '[bara]' '[°C]' '[-]' '[Nml min^-1]' ...
        '[W]' '[-]' '[-]' '[-]' '[-]' '[s]'};
    
    Response = Variables(result+3);
    Unit = Units(result+3);
    
    switch module
        case 1
            modulename = 'M1';
        case 2
            modulename = 'M2';
    end
    
    %% contourf plot
    qn = linspace(nmin, nmax, N);
    tn = linspace(nmin, nmax, N);
    [Qn,Tn] = meshgrid(qn,tn);
    Pn = polyval([mp tp], p)*ones(size(Qn));
    
    ypred = reshape(predict(mdl.(modulename).(Response{1}), ...
        [Qn(:) Pn(:) Tn(:)]), size(Qn));

    q = linspace(Qmin, Qmax, N);
    t = linspace(Tmin, Tmax, N);
    [Q,T] = meshgrid(q,t);
    
    if plotting == 1
        figure1 = figure('PaperSize',[19.98404194812 27.67743169791]);
        axes1 = axes('Parent',figure1,'FontSize',12,'FontName','Arial');        
        box(axes1,'on');
        hold(axes1,'all');
        [C,h] = contourf(Q, T, ypred, Ncontour, 'ShowText','off', 'LineWidth', 0.1,...
            'DisplayName', strcat(Response{1}, Unit{1}), 'HandleVisibility', 'off');
        v = [-0.3 0.3];
        clabel(C,h,v)
        caxis([0 max(ypred(:))])
        colormap gray
        c = colorbar;
        ylabel(c, strcat(Response, {' '},Unit));

        xlabel(axes1, {'Q_{LOHC} / ml·min^{-1}'},'FontSize',14,'FontName','Arial');
        ylabel(axes1, {'\vartheta / °C'},'FontSize',14,'FontName','Arial');
    end
    
    Q = Q(:);
    T = T(:);
    ypred = ypred(:);
end