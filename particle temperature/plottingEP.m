%load('EvaluationParticle2.mat')

legendnameQ = {'Raw', 'smooth', 'cumtrapz'};
legendnameT = {'280 °C', '290 °C', '300 °C', '310 °C'};
%legendnameX = {'MAX', '90 %', '80 %', '70 %', '60 %'};
legendnameProd = {'H', 'S'};
legendnamePdTX = {'Prod', 'dT', 'HG'};

plotQ = 0;                                                                  % QH2 depending on time
plotProdH = 0;                                                              % prod depending on time for H
plotProdS = 1;                                                              % prod depending on time for S
plotProdmax = 0;                                                            % maximum productivity H vs S at 310 degC
plotdTH = 0;                                                                % dT depending on time for H
plotdTS = 1;                                                                % dT depending on time for S  
plotPdTHG = 0;                                                              % P, dT and HG depending on time
plotProdd = 0;                                                              % Prod over d
plotdTd = 0;                                                                % dT over d
plotProddTd = 0;                                                            % Prod and dT over d
%% QH2 depending on time
if plotQ==1
  for i = 1:nd
     for j = 1:nTheta
       if j == 1 && ~isempty(H.t{i,j})
         figure('units','normalized','outerposition',[0 0 1 1])
         hold on
         title(sprintf('\\bf %s bei %d degC', strcat(Cat{1}, dstr{i}), Thetar(j)))
         h1 = plot(H.t{i,j}, H.QH2tRaw{i,j});
         h2 = plot(H.t{i,j}, H.QH2t{i,j});
         h3 = plot(H.t{i,j}, H.QH2{i,j});
         legend([h1 h2 h3], legendnameQ)
         hold off
       end  
       if ~isempty(S.t{i,j})
          if j == 1
            figure('units','normalized','outerposition',[0 0 1 1])
          end
          subplot(2,2,nTheta-j+1)
          hold on
          title(sprintf('\\bf %s bei %d degC', strcat(Cat{2}, dstr{i}), Thetar(j)))
  %        plot(S.t{i,j}, S.QH2tRaw{i,j},S.t{i,j}, S.QH2t{i,j},S.t{i,j}, S.QH2{i,j})
          s1 = plot(S.t{i,j}, S.QH2tRaw{i,j});     
          s2 = plot(S.t{i,j}, S.QH2t{i,j});
          s3 = plot(S.t{i,j}, S.QH2{i,j});
          A(j) = max(ylim);
          if j == 1
            legend([s1 s2 s3], legendnameQ)
          end
          if j >= 2
            if i == 6
              ymax = A(2);
            else 
              ymax = A(1);
            end
            ylim([0,ymax]);
          end
          hold off
       end
     end
  end
end

%% prod depending on time
if plotProdH == 1;
  j == 1;
  for i = 1:nd
    if i == 1;
      figure('units','normalized','outerposition',[0 0 1 1])
    end
    if ~isempty(H.t{i,j})
      hold on
      hh(i) = plot(H.t{i,j}, H.Prod{i,j});
      if i == nd
        title(sprintf('\\bf %s für %d degC', Cat{1}, Thetar(j)))
        legend(hh(hh~=0), dstr(hh~=0))
        xlim([0,360]);
        set(gca,'XGrid','on')
        set(gca,'YGrid','on')
        set(gca,'XTick',0:60:360);
        hold off
      end
    end
  end
end  

C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};                               % Cell array of colros.
if plotProdS == 1;  
  for j = 1:nTheta
     ss = zeros(1, nd);
     if j == 1;
         figure('units','normalized','outerposition',[0 0 1 1])
     end
     for i = 1:nd
       if ~isempty(S.t{i,j})
         subplot(2,2,nTheta-j+1)
         hold on
         title(sprintf('\\bf %s für %d degC', Cat{2}, Thetar(j)))
         ss(i) = plot(S.t{i,j}, S.Prod{i,j}, 'color', C{i});
         xlim([0,360]);
         set(gca,'XGrid','on')
         set(gca,'YGrid','on')
         set(gca,'XTick',0:60:360);
         
         A(i,j) = max(ylim);
         if i == nd
            ymax = max(A(i,1));
            ylim([0,ymax]);
         end
         hold off
         if i == nd && (j == 1 || j == 4)
           legend(ss(ss~=0), dstr(ss~=0))
         end       
      end
   end
 end
end


%% maximum productivity H vs S at 310 degC
j = 1;
ProdH = ones(nd,1)*-1;
ProdS = ones(nd,1)*-1;
if plotProdmax == 1
  for i = 1:nd
     if ~isempty(Results.H.Prod{i,j})
        ProdH(i) = Results.H.Prod{i,j}(1);
     end   
     ProdS(i) = Results.S.Prod{i,j}(1);
     
     if i == nd
       figure('units','normalized','outerposition',[0 0 1 1])
       hold on
       title(sprintf('maximum productivity at %d degC', Thetar(j)))
       h1 = plot(d(ProdH>0)*10^6, ProdH(ProdH>0), 'x');     
       s1 = plot(d*10^6, ProdS, 'x'); 
       
       ymax = max(ylim);
       legend([h1 s1], legendnameProd)
       ylim([0,ymax]);
       hold off
     end
  end
end  

%% dT depending on time
if plotdTH == 1;
  j == 1;
  for i = 1:nd
    if i == 1;
      figure('units','normalized','outerposition',[0 0 1 1])
    end
    if ~isempty(H.t{i,j})
      hold on
      hhh(i) = plot(H.t{i,j}, H.dT{i,j});
      if i == nd
        title(sprintf('\\bf %s für %d degC', Cat{1}, Thetar(j)))
        legend(hhh(hhh~=0), dstr(hhh~=0))
        xlim([0,360]);
        set(gca,'XGrid','on')
        set(gca,'YGrid','on')
        set(gca,'XTick',0:60:360);
        hold off
      end
    end
  end
end  

C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};                               % Cell array of colros.
if plotdTS == 1;  
  for j = 1:nTheta
     sss = zeros(1, nd);
     if j == 1;
         figure('units','normalized','outerposition',[0 0 1 1])
     end
     for i = 1:nd
       if ~isempty(S.t{i,j})
         subplot(2,2,nTheta-j+1)
         hold on
         title(sprintf('\\bf %s für %d degC', Cat{2}, Thetar(j)))
         sss(i) = plot(S.t{i,j}, S.dT{i,j}, 'color', C{i});
         xlim([0,360]);
         ylim([-2 2]);
         set(gca,'XGrid','on')
         set(gca,'YGrid','on')
         set(gca,'XTick',0:60:360);
         
         hold off
         if i == nd && (j == 1 || j == 4)
           legend(sss(sss~=0), dstr(sss~=0))
         end       
      end
   end
 end
end

%% % P, dT and X depending on time
if plotPdTHG == 1                                                               
  for i = 1:nd
     for j = 1:nTheta
       if j == 1 && ~isempty(H.t{i,j})
         figure('units','normalized','outerposition',[0 0 1 1])
         hold on
         title(sprintf('\\bf %s bei %d degC', strcat(Cat{1}, dstr{i}), Thetar(j)))
         [ax, h1, h3]= plotyy(H.t{i,j}, H.Prod{i,j},H.t{i,j}, H.HG{i,j}*100);
         lcolorh2 = get (gca, "ColorOrder");
         h2 = plot(H.t{i,j}, H.dT{i,j}, 'color', lcolorh2);
         set ([h1], "linestyle", "--");
         legend([h1 h2 h3], legendnamePdTX)
         ylabel (ax(1), "Prod and dT");
         ylabel (ax(2), "HG");
         
         xlim([0,360]);
         ylim([-2 2]);
         set(gca,'XGrid','on')
         set(gca,'YGrid','on')
         set(gca,'XTick',0:60:360);
         hold off
       end  
       if ~isempty(S.t{i,j})
          if j == 1
            figure('units','normalized','outerposition',[0 0 1 1])
          end
          subplot(2,2,nTheta-j+1)
          hold on
          title(sprintf('\\bf %s bei %d degC', strcat(Cat{2}, dstr{i}), Thetar(j)))
          [ax, s1, s3]= plotyy(S.t{i,j}, S.Prod{i,j},S.t{i,j}, S.HG{i,j}*100);
          lcolors2 = get (gca, "ColorOrder");
          s2 = plot(S.t{i,j}, S.dT{i,j}, 'color', lcolors2);
          set ([s1], "linestyle", "--");
          if j == 1
            legend([s1 s2 s3], legendnamePdTX)
            ylabel (ax(1), "Prod and dT");
            ylabel (ax(2), "HG");
          end
          xlim([0,360]);
          ylim([-2 2]);
          set(gca,'XGrid','on')
          set(gca,'YGrid','on')
          set(gca,'XTick',0:60:360);
          hold off
       end
     end
  end
end

%% prod over d
if plotProdd == 1 || plotProddTd == 1
  ProdS = ones(nd, nTheta)*-1;
  for j = 1:nTheta
    for i = 1:nd
       if ~isempty(Results.S.Prod{i,j})
          ProdS(i,j) = Results.S.Prod{i,j}(1);
       end   
       
%       if i == nd
%         figure('units','normalized','outerposition',[0 0 1 1])
%         hold on
%         title(sprintf('maximum productivity at %d degC', Thetar(j)))
%         h1 = plot(d(ProdH>0)*10^6, ProdH(ProdH>0), 'x');     
%         s1 = plot(d*10^6, ProdS, 'x'); 
%         
%         ymax = max(ylim);
%         legend([h1 s1], legendnameProd)
%         ylim([0,ymax]);
%         hold off
%       end
    end
  end  
  if plotProdd == 1 
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    title(sprintf('prod over d'))
    s1 = plot(d(ProdS(:,1)>0)*10^6, ProdS(ProdS(:,1)>0,1), '--o');
    s2 = plot(d(ProdS(:,2)>0)*10^6, ProdS(ProdS(:,2)>0,2), '--o');
    s3 = plot(d(ProdS(:,3)>0)*10^6, ProdS(ProdS(:,3)>0,3), '--o');
    s4 = plot(d(ProdS(:,4)>0)*10^6, ProdS(ProdS(:,4)>0,4), '--o');
    legend([s1 s2 s3 s4], legendnameT)
  end
end  

%% dT over d
if plotdTd == 1 || plotProddTd == 1
  dTS = ones(nd, nTheta)*-5;
  for j = 1:nTheta
    for i = 1:nd
       if ~isempty(Results.S.dT{i,j})
          dTS(i,j) = Results.S.dT{i,j}(1);
       end   
       
%       if i == nd
%         figure('units','normalized','outerposition',[0 0 1 1])
%         hold on
%         title(sprintf('maximum productivity at %d degC', Thetar(j)))
%         h1 = plot(d(ProdH>0)*10^6, ProdH(ProdH>0), 'x');     
%         s1 = plot(d*10^6, ProdS, 'x'); 
%         
%         ymax = max(ylim);
%         legend([h1 s1], legendnameProd)
%         ylim([0,ymax]);
%         hold off
%       end
    end
  end  
  if plotdTd == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    title(sprintf('dT over d'))
    s5 = plot(d(dTS(:,1)>-5)*10^6, dTS(dTS(:,1)>-5,1), '-x');
    s6 = plot(d(dTS(:,2)>-5)*10^6, dTS(dTS(:,2)>-5,2), '-x');
    s7 = plot(d(dTS(:,3)>-5)*10^6, dTS(dTS(:,3)>-5,3), '-x');
    s8 = plot(d(dTS(:,4)>-5)*10^6, dTS(dTS(:,4)>-5,4), '-x');
    legend([s5 s6 s7 s8], legendnameT)
  end  
end  

%% Prod and dT over d
if plotProddTd == 1
%  figure('units','normalized','outerposition',[0 0 1 1])
%  hold on
%  title(sprintf('Prod and dT over d'))
%  s1 = plot(d(ProdS(:,1)>0)*10^6, ProdS(ProdS(:,1)>0,1), '--o');
%  s2 = plot(d(ProdS(:,2)>0)*10^6, ProdS(ProdS(:,2)>0,2), '--o');
%  s3 = plot(d(ProdS(:,3)>0)*10^6, ProdS(ProdS(:,3)>0,3), '--o');
%  s4 = plot(d(ProdS(:,4)>0)*10^6, ProdS(ProdS(:,4)>0,4), '--o');
%  set(gca,'ColorOrderIndex',1);
%  s5 = plot(d(dTS(:,1)>0)*10^6, dTS(dTS(:,1)>0,1), '-x');
%  s6 = plot(d(dTS(:,2)>0)*10^6, dTS(dTS(:,2)>0,2), '-x');
%  s7 = plot(d(dTS(:,3)>0)*10^6, dTS(dTS(:,3)>0,3), '-x');
%  s8 = plot(d(dTS(:,4)>0)*10^6, dTS(dTS(:,4)>0,4), '-x');
%  legend([s1 s2 s3 s4], legendnameT)
  
   figure('units','normalized','outerposition',[0 0 1 1])
   hold on
   title(sprintf('Prod and dT over d'))
   [ax, s1, s5] = plotyy(d(ProdS(:,1)>0)*10^6, ProdS(ProdS(:,1)>0,1), d(dTS(:,1)>-5)*10^6, dTS(dTS(:,1)>-5,1));
   [~, s2, s6] = plotyy(d(ProdS(:,2)>0)*10^6, ProdS(ProdS(:,2)>0,2), d(dTS(:,2)>-5)*10^6, dTS(dTS(:,2)>-5,2));
   [~, s3, s7] = plotyy(d(ProdS(:,3)>0)*10^6, ProdS(ProdS(:,3)>0,3), d(dTS(:,3)>-5)*10^6, dTS(dTS(:,3)>-5,3));
   [~, s4, s8] = plotyy(d(ProdS(:,4)>0)*10^6, ProdS(ProdS(:,4)>0,4), d(dTS(:,4)>-5)*10^6, dTS(dTS(:,4)>-5,4));
   set ([s1 s2 s3 s4], 'linestyle', '--', 'marker', 'o');
   set ([s5 s6 s7 s8], 'linestyle', '-', 'marker', 'x');
   set ([s1 s5], 'color', C{1})
   set ([s2 s6], 'color', C{2})
   set ([s3 s7], 'color', C{3})
   set ([s4 s8], 'color', C{4})
   set(ax,{'ycolor'},{'k';'k'})
   xlim([0,200]);
   legend([s1 s2 s3 s4], legendnameT)
   ylabel (ax(1), 'Prod');
   ylabel (ax(2), 'dT');
   set(ax(2),'YLim',[-2 2])
end 