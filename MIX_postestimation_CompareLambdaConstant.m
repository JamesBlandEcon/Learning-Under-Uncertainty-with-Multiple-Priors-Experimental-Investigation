clear; clc;
set(0,'defaulttextInterpreter','latex') %latex axis labels

%% 
load 'C:\Users\jbland\Documents\LearningSims\MixBetaLambdaConstant'
  


% Remove BurnIn
varlist = {'THETA','THETA_MU','THETA_V','PRTYPE','ITYPE','MIXPROBS','LAMBDA'};
for vv = 1:numel(varlist)
    str = [varlist{vv} '(:,:,1:BurnIn)=[];'];
    eval(str)
end


% save some of the simuation

LC.lambda = exp(LAMBDA);
LC.gamma  = THETA(:,1,:);
LC.gammaM = THETA_MU(1,1,:);
LC.gammaV = THETA_V(1,1,:);

clearvars -except LC

%%

load 'C:\Users\jbland\Documents\LearningSims\MixBeta'

% Remove BurnIn
varlist = {'THETA','THETA_MU','THETA_V','PRTYPE','ITYPE','MIXPROBS'};
for vv = 1:numel(varlist)
    str = [varlist{vv} '(:,:,1:BurnIn)=[];'];
    eval(str)
end


% save some of the simuation

LH.lambda = exp(THETA(:,2,:));
LH.gamma  = THETA(:,1,:);
LH.gammaM = THETA_MU(1,1,:);
LH.gammaV = THETA_V(1,1,:);
LH.lambdaM = THETA_MU(1,2,:);
LH.lambdaV = THETA_V(2,2,:);

clearvars -except LC LH

%%

lw = 1.5;
rng(42)
ci = 0.9
hsf = figure('units','inches','Position',[0.1 1 8 8]);
FFXX  = linspace(0,1,101);
h = subplot(2,2,1)
    hold all
        data = LH.lambda;
        plotscript
        delete(l4)
        xlims = [0 max(xm)];
        xlim(xlims)
        FFXX = linspace(xlims(1),xlims(2),1001);
        
        cdfL = prctile(permute(logncdf(FFXX,LH.lambdaM,sqrt(LH.lambdaV)),[2 3 1]),100*[(1-ci)/2 1-(1-ci)/2],2);
        pltcdf = plot(FFXX,cdfL,'-b','LineWidth',lw);
        
        %cdfLC = prctile(permute(logncdf(FFXX,LC.lambdaM,sqrt(LC.lambdaV)),[2 3 1]),100*[(1-ci)/2 1-(1-ci)/2],2);
        % pltcdf = plot(FFXX,cdfLC,'-r','LineWidth',lw);
        pltLambdaCons = plot([1 1]*prctile(LC.lambda,50),[0 1],'--r','LineWidth',lw);
                        %plot([1 1]*prctile(LC.lambda,100-100*ci/2),[0 1],'-r','LineWidth',2)
        
        %legend([l1 l2 pltcdf(1)],{'Ind 50th percentile','Ind 5th-95th percentile','Pop 5-95th percentile'},'Interpreter','latex','Location','best','Box','off')               
        ylabel('Cumulative probability')       
        xlabel('$\lambda$')
        title('(a) Choice precision')
    hold off
h = subplot(2,2,2); % compare cdfs of gamma
    hold all
        gmH = median(LH.gamma,3);
        gmHsorted = sort(gmH);
        gmC = median(LC.gamma,3);
        gmCsorted = sort(gmC);
        xlims = [min(gmHsorted) max(gmHsorted)];
        xlim(xlims)
        FFXX = linspace(xlims(1),xlims(2),1001);
        cdfGH  = permute(prctile(normcdf(FFXX,LH.gammaM,sqrt(LH.gammaV)),100*[(1-ci)/2 1-(1-ci)/2],3),[3 2 1]);
        pltcdfH = plot(FFXX,cdfGH,'-b','LineWidth',lw);
        cdfGC  = permute(prctile(normcdf(FFXX,LC.gammaM,sqrt(LC.gammaV)),100*[(1-ci)/2 1-(1-ci)/2],3),[3 2 1]);
        pltcdfC = plot(FFXX,cdfGC,'--r','LineWidth',lw);
        
        plot(gmHsorted,linspace(0,1,numel(gmHsorted)),'.b')
        plot(gmCsorted,linspace(0,1,numel(gmHsorted)),'xr')
        title('(b) Risk aversion')
        xlabel('$\gamma$')
    hold off
h = subplot(2,2,3); % scatter plot of risk aversion
    hold all
        plot(gmH,gmC,'.k')
        plot(gmHsorted,gmHsorted,':k')
        xlim([-1 1])
        ylim([-1 1])
        xlabel('$\gamma$ assuming heterogeneous $\lambda$')
        ylabel('$\gamma$ assuming homogeneous $\lambda$')
        title('(c) Sensitivity to $\lambda$')
    hold off
h = subplot(2,2,4); %legend
    plt.XTick = [];
    plt.YTick = [];
    plt.XColor = [1 1 1];
    plt.YColor = [1 1 1];
    hold all
    axis off
        p1 = plot(-1,-1,'.k')
        p1a = plot(-1,-1,'xr')
        p2 = plot(-1,-1,'-k')
        p3 = plot(-1,-1,'-b','LineWidth',lw)
        p4 = plot(-1,-1,'--r','LineWidth',lw)
        p5 = plot(-1,-1,':k')
        ylim([0 1])
        xlim([0 1])
        leg = legend([p1 p1a p2 p3 p4 p5],...
        {[' ' newline 'Individual 50th percentile', newline,'(heterogeneous $\lambda$)' newline ' '], ...
        [' ' newline 'Individual 50th percentile', newline,'(homogeneous $\lambda$)' newline ' '], ...
        [' ' newline 'Individual 5th-95th percentile' newline ' '],...
         [' ' newline 'Population 5th-95th percentile' newline '(heterogeneous $\lambda$)' newline ' '],...
         [' ' newline 'Population 5th-95th percentile' newline '(homogeneous $\lambda$)' newline ' '],... 
         [' ' newline '$45^\circ$ line' newline ' ']},...
         'Box','off','Interpreter','latex','FontSize',11)
saveas(hsf,'figures/heterogeneousLambdaCompare.png')
%     hold all
%         dll = log(median(LH.lambda,3))-log(median(LC.lambda));
%         plot(dll,gmH-gmC,'.k')
%         ylim([-1 1])
%         xlim([min(dll) max(dll)])
%         xlabel('$\log\lambda_i-\log\bar\lambda$')
%         ylabel('$\gamma_i-\bar\gamma_i$')
%         title('(d) Sensitivity to $\lambda$')
%     hold off

%% Fraction risk averse

RiskAverseH = 1-normcdf(0,LH.gammaM(:),sqrt(LH.gammaV(:)));
RiskAverseC = 1-normcdf(0,LC.gammaM(:),sqrt(LC.gammaV(:)));

h = figure;
hold all
    [F,XI]=ksdensity(RiskAverseH);
    fH = plot(XI,F,'-b')
    [F,XI]=ksdensity(RiskAverseC);
    fC = plot(XI,F,'-r')
    legend([fH fC],{'Heterogeneous $\lambda$','Homogeneous $\lambda$'},'Location','best','box','off','Interpreter','latex')
    xlabel('Fraction risk averse')
    ylabel('Posterior density')
hold off

disp(['Hetrogeneous: ' num2str(mean(RiskAverseH)) ' (' num2str(std(RiskAverseH)) ')'])
disp(['Homogeneous:  ' num2str(mean(RiskAverseC)) ' (' num2str(std(RiskAverseC)) ')'])










