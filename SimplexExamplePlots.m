clear; clc;

FontSize = 20;
set(0,'defaultLegendFontSize',FontSize)
set(0,'defaultAxesFontSize',FontSize)
set(0,'defaultLegendInterpreter','latex')


delta = 0.1;
alpha = 0.4;
S = [0 0 1 1 1 1];
%S = zeros(1,100);

nic = 6; dic = 0.15;
SS = [0 3 6];
for kk = 1:numel(SS)
    
    ss = SS(kk);
    h = figure
    set(gca,'TickLabelInterpreter','none')
    
    hold all
    
    urnB = [1 3 2];

        
        % boundary of Simplex 
        l1 = patch([0 1 0 0],[1 0 0 1],'white','LineWidth',2,'EdgeColor','k','LineStyle','--');
        
        
        
        
        % boundary of set of priors
        patchThis = [[delta 1-2.*delta delta delta]' [delta delta 1-2*delta delta]'];
        l2 = patch(patchThis(:,1),patchThis(:,2),'white','LineWidth',2,'EdgeColor','b');
        
        %area([delta 1-2*delta],[delta delta], 'EdgeColor',[1 1 1],'FaceColor',[1,1,1]);
        
        

        
        %fplot(@(x) (x-delta) + 1-2*delta,[0 1],'--g')
        %fplot(@(x) (x-(1-2.*delta)) + delta,[0 1],'--g')
        
         
            b = sum(S(1:ss));
            w = sum(1-S(1:ss));
            
        
        xlim([0 1])
        ylim([0 1])
            
        if ss>0
            L = zeros(1,3);
            for ll = 1:3
                L(ll) = ll.^b.*(4-ll).^w.*(1-3.*delta)+delta.*sum((1:3).^b.*(4-1:3).^w);
            end

            [Lmax,II] = max(L);
            if II == 1
                mlp = [1-2*delta, delta];
            elseif II==2
                mlp = [delta, delta];
            else
                mlp = [delta, 1-2*delta];
            end
            [Lmax II mlp]
            l4 = plot(mlp(1),mlp(2),'ok','MarkerFaceColor','k','MarkerSize',10);
            intercept = (alpha.*Lmax-2.^(w+b))./(3.^b-2.^(w+b));
            slope = -(3.^w-2.^(w+b))./(3.^b-2.^(w+b));
            xx = linspace(0,1,1001);
            yy = intercept+slope.*xx;
            
            l5 = plot(xx((xx+yy)<=1),yy((xx+yy)<=1),'-.r','LineWidth',3 ); %fplot(@(x) intercept + slope.*x,[0 1],'-r','LineWidth',2);
            if II == 1
                patchThis = [delta delta;...
                            mlp;...
                            1-delta-(intercept+slope.*(1-delta))./(1+slope) (intercept+slope.*(1-delta))./(1+slope);...
                            delta intercept+slope.*delta;...
                            delta delta];
            elseif II == 3
                patchThis = [delta delta;...
                            mlp;...
                            1-delta-(intercept+slope.*(1-delta))./(1+slope) (intercept+slope.*(1-delta))./(1+slope);...
                            (delta-intercept)./slope delta];    
            end
            l6 = patch(patchThis(:,1),patchThis(:,2),'b','FaceAlpha',0.2,'LineStyle','none');
            
            % set of posterior beliefs
            PriorSet = [patchThis 1-sum(patchThis,2)];
            PostSet  = NaN.*zeros(size(PriorSet));
            
            
            PostSet = repmat((urnB./4).^b.*(1-urnB./4).^w,[size(PostSet,1) 1]).*PriorSet;
            PostSet = PostSet./repmat(sum(PostSet,2),[1 3]);
            l7 = patch(PostSet(:,1),PostSet(:,2),'r','FaceAlpha',0.6,'LineStyle','-','EdgeColor','r');
            PrBlack = sum(PostSet.*repmat(urnB./4,[size(PostSet,1) 1]),2);
            [WCS,WCSii] = min(PrBlack);
            l8 = plot(PostSet(WCSii,1),PostSet(WCSii,2),'xk','MarkerSize',10,'LineWidth',3);
            
            PriorBlack = sum(PriorSet.*repmat(urnB./4,[size(PriorSet,1) 1]),2);
        else
            %l4 = [];
            l5 = [];
            l6 = [];
            l7 = [];
            l8 = [];
            PriorSet = [patchThis 1-sum(patchThis,2)];
            PostSet  = NaN.*zeros(size(PriorSet));
            PostSet = repmat((urnB./4).^b.*(1-urnB./4).^w,[size(PostSet,1) 1]).*PriorSet;
            PostSet = PostSet./repmat(sum(PostSet,2),[1 3]);
            PrBlack = sum(PostSet.*repmat(urnB./4,[size(PostSet,1) 1]),2);
            
            PriorBlack = sum(PriorSet.*repmat(urnB./4,[size(PriorSet,1) 1]),2);
            l4 = patch([delta 1-2.*delta delta delta]', [delta delta 1-2*delta delta]','b','FaceAlpha',0.3,'LineStyle','none');
            [WCS,WCSii] = min(PrBlack);
            l8 = plot(PostSet(WCSii,1),PostSet(WCSii,2),'xk','MarkerSize',10,'LineWidth',3);
        end
        
        % indiffence curves
        xx = linspace(0,0.5,101);
        l3 = plot(xx ,xx,'--k');
        for ii = 1:nic
            fplot(@(x) x + ii.*dic,[0 (1-ii.*dic)./2],'--k');
            fplot(@(x) x - ii.*dic,[0 (1+ii.*dic)./2],'--k');
        end
        
        %title(['w = ' int2str(w) ', b = ' int2str(b)])
        xl = xlabel('$\mu_1$','Interpreter','latex')
        yl = ylabel('$\mu_3$','Interpreter','latex')
        %xl.FontSize = FontSize;
        %yl.FontSize = FontSize;
        xlim([0 1])
        ylim([0 1])
        if kk==1
            % get legend patch opacity correct
            leg = {'Simplex','Indifference curves ~', 'Initial priors ($M_0$)','Likely priors ($M_1$)'};
            [l,icons,plots,txt] = legend([l1 l3 l2 l4],leg)
            l.Box='off';
            PatchInLegend = findobj(icons, 'type', 'patch');
            set(PatchInLegend, 'facea', 0.5)
            for ii = 1:4
                icons(ii).Interpreter = 'latex';
                icons(ii).FontSize = FontSize;
            end
            
        end
        if kk==2
            % get legend patch opacity correct
            [l,icons,plots,txt] = legend([ l4 l5 l7 l8],'Most likely prior in $M_0$','$L = \alpha L^*$','Set of posterior beliefs','Worst-case scenario')
            l.Box='off';
            PatchInLegend = findobj(icons, 'type', 'patch');
            set(PatchInLegend, 'facea', 0.5)
            for ii = 1:4
                icons(ii).Interpreter = 'latex';
                icons(ii).FontSize = FontSize;
            end
        end
        if kk==3
            % do nothing
        end
        if kk>1
            counter = 1;
            dotoffset = 0.025;
        for ll = (SS(1)+1):SS(kk) %for ll = (SS(kk-1)+1):SS(kk)
            if S(ll) == 1;
                dot = 'k';
            else
                dot = 'w';
            end
            plot(0.3-dotoffset+dotoffset.*(counter-2),0.95,'ok','MarkerFaceColor',dot)
            counter = counter+1;
        end
        end

        

        %set(PatchInLegend,'f')
        % remove whitespace around figure
         ax = gca;
         outerpos = ax.OuterPosition;
         ti = ax.TightInset; 
         left = outerpos(1) + ti(1);
         bottom = outerpos(2) + ti(2);
         ax_width = outerpos(3) - ti(1) - ti(3);
         ax_height = outerpos(4) - ti(2) - ti(4);
         ax.Position = [left bottom ax_width ax_height];
         
         
    hold off
    saveas(h,['figures/SimplexExample_S' int2str(ss) '.png'])

    PrBlack_Post{kk} = PrBlack;
    
    PrBlack_Prior{kk} = PriorBlack;
    if kk>1;
        mlpBlack{kk} = sum(urnB./4.*[mlp (1-sum(mlp))]);
    end
    
end

%%

close all
FontSize = 16;
set(0,'defaultLegendFontSize',FontSize)
set(0,'defaultAxesFontSize',FontSize)
set(0,'defaultLegendInterpreter','latex')



h = figure;
offset = 0.3;
ha = axes;
hold all

ha.XColor = [1 1 1];
xl = xlabel('Draw');
xl.Color = [0 0 0];
xl.FontSize = FontSize;


for kk = 1:numel(SS)
        
        b = sum(S(1:SS(kk)));
        w = sum(1-S(1:SS(kk)));
        % Set of priors
        l1 = plot(SS(kk).*[1 1]-offset,[min(PrBlack_Prior{kk}) max(PrBlack_Prior{kk})],'--b','LineWidth',2)
        % Set of posterior beliefs
        l2 = plot(SS(kk).*[1 1]+offset,[min(PrBlack_Post{kk}) max(PrBlack_Post{kk})],'-r','LineWidth',2)
        % worst-case scenatior
        l4 = plot(SS(kk)+offset,min(PrBlack_Post{kk}),'xk','MarkerSize',10,'LineWidth',3);
        if kk>1
            % most likely prior
            l3 = plot(SS(kk)-offset,mlpBlack{kk},'ok','MarkerFaceColor','k','MarkerSize',10)
        end
        if kk>1
            l6 = plot(SS(kk),b./(b+w),'+k','MarkerSize',10,'LineWidth',3);
        end
        
        if kk>1
            counter = 1;
            dotoffset = 0.2;
        for ll = (SS(kk-1)+1):SS(kk)
            if S(ll) == 1;
                dot = 'k';
            else
                dot = 'w';
            end
            plot(SS(kk)+dotoffset.*(counter-2),0.95,'ok','MarkerFaceColor',dot)
            counter = counter+1;
        end
        end
end
    % allowable range
    l5 = plot([-1 max(SS)+1],[0.75 0.75],'-.k','LineWidth',2)
    plot([-1 max(SS)+1],[0.25 0.25],'-.k','LineWidth',2)
    % indifference curves
    for yy = 0.1:0.1:0.9;
        l7 = plot([-1 max(SS)+1],yy.*[1 1],'--k')
    end
    legend([l1 l2 l3 l5 ],'Likely priors ($M_1$)','Set of posterior beliefs','Most likely prior in $M_0$','Allowable range','Location','SouthWest','Interpreter','latex')
    ylim([0 1])
    xlim([-1 max(SS)+1])
    yl = ylabel('Pr[black]')
    yl.FontSize = FontSize;
    
    
    set(ha,'XTick',[]);
    
    %set(ha,'XTick',SS);
    % remove whitespace around figure
         ax = gca;
         outerpos = ax.OuterPosition;
         ti = ax.TightInset; 
         left = outerpos(1) + ti(1);
         bottom = outerpos(2) + ti(2);
         ax_width = outerpos(3) - ti(1) - ti(3);
         ax_height = outerpos(4) - ti(2) - ti(4);
         ax.Position = [left bottom ax_width ax_height];
         
hold off
saveas(h,['figures/SimplexExample_PrBlack.png'])

%% Beta priors

N0 = 10;
% use same range of p0 as in Simplex priors
p0L =min(PrBlack_Prior{1});
p0R =max(PrBlack_Prior{1});

% parameters from Yaroslav's original Beta priors example
% p0L = 0.45;
% p0R = 0.55;
% N0 = 4;
% alpha = 0.95

h = figure;
ha = axes;
hold all
    ha.XColor = [1 1 1];
    xl = xlabel('Draw');
    xl.Color = [0 0 0];
    xl.FontSize = FontSize;
    
    prior = [p0L p0R];
    post = [p0L p0R];
    Lprior = NaN.*[1 1];
    for kk = 1:numel(SS)
        if kk>1
            ss = SS(kk);
            b = sum(S(1:ss))
            w = ss-b
            % Find mlp
            pgrid = linspace(p0L,p0R,1001);
            LikePrior = beta(b+N0.*pgrid,w+N0.*(1-pgrid))./beta(N0.*pgrid,N0.*(1-pgrid));
            [Lmax,II] = max(LikePrior);
            mlp = pgrid(II);
            %mlp = sum((Lprior==max(Lprior)).*prior);
            % most likly prior
            
           
            
            %pAlpha = fsolve(L,0.5);
            prior = pgrid(LikePrior>=(alpha.*Lmax));
            prior = [min(prior) max(prior)];
            for pp = 1:numel(prior)
                post(pp) = (prior(pp).*N0+b)./(N0+b+w)
            end
            l6 = plot(SS(kk),b./(b+w),'+k','MarkerSize',10,'LineWidth',3);
            
        end
        % set of priors
        l1 = plot((SS(kk)-offset).*[1 1],prior,'--b','LineWidth',2)
        % set of posteriors
        l2 = plot((SS(kk)+offset).*[1 1],post,'-r','LineWidth',2)
        l4 = plot((SS(kk)+offset),min(post),'xk','MarkerSize',10,'LineWidth',3);
        
        
            
        
        if kk>1
            l3 =  plot((SS(kk)-offset),mlp,'ok','MarkerFaceColor','k','MarkerSize',10)
            counter = 1;
            dotoffset = 0.2;
        for ll = (SS(kk-1)+1):SS(kk)
            if S(ll) == 1;
                dot = 'k';
            else
                dot = 'w';
            end
            plot(SS(kk)+dotoffset.*(counter-2),0.95,'ok','MarkerFaceColor',dot)
            counter = counter+1;
        end
        end
        
    end
    for yy = 0.1:0.1:0.9;
        l7 = plot([-1 max(SS)+1],yy.*[1 1],'--k')
    end
    l5 = plot([-1 max(SS)+1],[1 1],'-.k','LineWidth',2)
    plot([-1 max(SS)+1],[0 0],'-.k','LineWidth',2)
    legend([l4 l6 l7 ],'Worst-case scenario','b/(b+w)','Indifference curves','Location','SouthEast')
    ylim([0 1])
    xlim([-1 max(SS)+1])
    yl = ylabel('Pr[black]')
    yl.FontSize = FontSize;
    
    set(ha,'XTick',[]);
    
    
    % The Beta priors diagram covers more space, so remove legend here if
    % we only want to show one
    %legend([l1 l2 l3 l4 l5 l6],'Likely priors (M_1)','Set of posterior beliefs','Most likely prior in M_0','Worst-case scenario','Allowable range of Pr[black]','b/(b+w)','Location','SouthEast')
    % remove whitespace around figure
         ax = gca;
         outerpos = ax.OuterPosition;
         ti = ax.TightInset; 
         left = outerpos(1) + ti(1);
         bottom = outerpos(2) + ti(2);
         ax_width = outerpos(3) - ti(1) - ti(3);
         ax_height = outerpos(4) - ti(2) - ti(4);
         ax.Position = [left bottom ax_width ax_height];
hold off
saveas(h,['figures/BetaExample_PrBlack.png'])
%saveas(h,['BetaExample_PrBlack_N0_' int2str(N0) '.png'])


