clear; clc;


%load 'C:\Users\jbland\Documents\LearningSims\MixBeta'
%load 'C:\Users\jbland\Documents\LearningSims\MixBetaRestricted'
%load 'C:\Users\jbland\Documents\LearningSims\MixBetaExpt2'
load 'C:\Users\jbland\Documents\LearningSims\MixBeta_Both' 

set(0,'defaulttextInterpreter','latex') %latex axis labels
%%

TypeList = {'A-B C-B','A-B C-MP','A-MP C-B','A-MP C-MP'};


if     strcmp(modelname,'MixBeta') || strcmp(modelname,'MixBeta_Both') || strcmp(modelname,'MixBetaExpt2NoSymmetric') || strcmp(modelname,'MixBetaRestricted')
    ParamLabels = {'\gamma','log\lambda','A-B-N_0','A-B-p_0','A-MP-\alpha','A-MP-N_0','A-MP-p_{0mid}','A-MP-p_{0spread}','C-B-N_0','C-B-p_0','C-MP-\alpha','C-MP-N_0','C-MP-p_{0mid}','C-MP-p_{0spread}'};
    ParamTxt =    {'gamma','lambda'     ,'A-B-N0' ,'A-B-p0','A-MP-alpha','A-MP-N0','A-MP-p0mid','A-MP-p0spread','C-B-N0','C-B-p0','C-MP-alpha','C-MP-N0','C-MP-p0mid','C-MP-p0spread'};
    transform =   {''     ,''           ,'exp',   'normcdf','normcdf'   ,'exp'    ,'normcdf'   ,'normcdf'      ,'exp',   'normcdf','normcdf'   ,'exp'    ,'normcdf'   ,'normcdf'      ,                 };
    PlotThis1 = [5 11; 6 12; 7 13; 8 14];
    PlotThisB = [3 4 9 10];
    PlotThisMP = [5 6 7 8 11 12 13 14];
    xlimsB    = [0 40; 0 1; 0 40; 0 1];
    xlimsMP    = [0 1; 0 40; 0 1; 0 1; 0 1; 0 40; 0 1; 0 1];
    PlotRange1 = [-5 log(30);  0.4 0.6; -2 2; -2 2];
    ParamTable = {'\gamma','\lambda','N_{0,AB}','p_{0,AB}','\alpha_{AM}','N_{0,AM}','p_{0\text{m},AM}','p_{0\text{s},AM}','N_{0,CB}','p_{0,CB}','\alpha_{CM}','N_{0,CM}','p_{0\text{m},CM}','p_{0\text{s},CM}'};
elseif strcmp(modelname,'MixBetaExpt2') || strcmp(modelname,'MixBetaExpt2short')
    ParamLabels = {'\gamma','log\lambda','A-B-N_0','A-B-p_0','A-MP-\alpha','A-MP-N_0','A-MP-p_{0L}','C-B-N_0','C-B-p_0','C-MP-\alpha','C-MP-N_0','C-MP-p_{0L}'};
    ParamTxt =    {'gamma','lambda'     ,'A-B-N0' ,'A-B-p0','A-MP-alpha','A-MP-N0','A-MP-p0L','C-B-N0','C-B-p0','C-MP-alpha','C-MP-N0','C-MP-p0L'};
    transform =   {''     ,''           ,'exp',   'normcdf','normcdf'   ,'exp'    ,'0.5.*normcdf' ,'exp',   'normcdf','normcdf'   ,'exp'    ,'0.5.*normcdf'   };
    PlotThis1 = [5 11; 6 12; 7 13; 8 14];
    PlotThisB = [3 4 9 10];
    PlotThisMP = [5 6 7 10 11 12];
    xlimsB    = [0 40; 0 1; 0 40; 0 1];
    xlimsMP    = [0 1; 0 40; 0 0.5; 0 1; 0 10; 0 0.5];
    PlotRange1 = [-5 log(30);  0.4 0.6; -2 2; -2 2];
    ParamTable = {'\gamma','\lambda','N_{0,AB}','p_{0,AB}','\alpha_{AM}','N_{0,AM}','p_{0\text{m},AM}','p_{0\text{s},AM}','N_{0,CB}','p_{0,CB}','\alpha_{CM}','N_{0,CM}','p_{0\text{m},CM}','p_{0\text{s},CM}'};
elseif strcmp(modelname,'MixSimplex')
    ParamLabels = {'\gamma','log\lambda','A-B-p_1','A-B-p_2','A-MP-\alpha','A-MP-p_{min}','C-B-p_1','C-B-p_2','C-MP-\alpha','C-MP-p_{min}'};
    ParamTxt =    {'gamma','lambda'     ,'A-B-p1' ,'A-B-p2','A-MP-alpha',  'A-MP-pmin',   'C-B-p1', 'C-B-p2', 'C-MP-alpha', 'C-MP-pmin'};
    transform =   {''     ,''           ,'normcdf','normcdf','normcdf'   ,'1./3.*normcdf','normcdf','normcdf','normcdf','1./3*normcdf'};
    PlotThis1 = [3 7; 4 8; 5 9; 6 10];
    PlotThisB = [3 4 7 8];
    PlotThisMP = [5 6 9 10];
    xlimsB    = [0 1; 0 1; 0 1; 0 1];
    xlimsMP    = [0 1; 0 1/3; 0 1; 0 1/3];
    PlotRange1 = [0 1;  0 1; 0 1; 0 1];
    ParamTable = {'\gamma','\lambda','p_{1AB}','p_{2AB}','\alpha_{AM}','p_{\text{min},AM}','p_{1CB}','p_{2CB}','\alpha_{CM}','p_{\text{min},CM}'};
elseif strcmp(modelname,'MixSimplexGeneral')
    ParamLabels = {'\gamma','log\lambda','A-B-p_1','A-B-p_2','A-MP-\alpha','A-MP-p_{min1}','A-MP-p_{min2}','A-MP-p_{min3}','C-B-p_1','C-B-p_2','C-MP-\alpha','C-MP-p_{min}','C-MP-p_{min2}','C-MP-p_{min3}'};
    ParamTxt =    {'gamma','lambda'     ,'A-B-p1' ,'A-B-p2','A-MP-alpha',  'A-MP-pmin1','A-MP-pmin2','A-MP-pmin3',   'C-B-p1', 'C-B-p2', 'C-MP-alpha', 'C-MP-pmin1','C-MP-pmin2','C-MP-pmin3'};
    transform =   {''     ,''           ,'normcdf','normcdf','normcdf','normcdf','normcdf','normcdf','normcdf','normcdf','normcdf','normcdf','normcdf','normcdf'};
    PlotThis1 = [3 7; 4 8; 5 9; 6 10];
    PlotThisB = [3 4 9 10];
    PlotThisMP = [5 6 7 8 11 12 13 14];
    xlimsB    = [0 1; 0 1; 0 1; 0 1];
    xlimsMP    = [0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1];
    PlotRange1 = [0 1;  0 1; 0 1; 0 1];
    ParamTable = {'\gamma','\lambda','p_{1AB}','p_{2AB}','\alpha_{AM}','p_{\text{min1},AM}','p_{\text{min2},AM}','p_{\text{min3},AM}','p_{1CB}','p_{2CB}','\alpha_{CM}','p_{\text{min1},CM}','p_{\text{min2},CM}','p_{\text{min3},CM}'};
end




%% Remove BurnIn
varlist = {'THETA','THETA_MU','THETA_V','PRTYPE','ITYPE','MIXPROBS'};

%BurnIn=5000; % override from initial simulation

for vv = 1:numel(varlist)
    str = [varlist{vv} '(:,:,1:BurnIn)=[];'];
    eval(str)
end


%% Diagnostics: have we got to a stationary distribution yet?
MP = permute(MIXPROBS,[3 2 1]);
h = figure;
for mm = 1:size(MP,2)
    subplot(size(MP,2),1,mm)
    hold all
    title(TypeList{mm})
    plot(MP(:,mm))
    hold off
end

TM = permute(THETA_MU,[3 2 1]) ;
h = figure;
for vv =1:size(TM,2)
    subplot(size(TM,2)/2,2,vv)
    hold all
    %title(TypeList{mm})
    plot(TM(:,vv))
    hold off
end



%% MixProbs posterior 

MP = permute(MIXPROBS,[3 2 1]);
xx = linspace(0,1,201);
h = figure;
hold all
for tt = 1:numel(TypeList)
    yy = ksdensity(MP(:,tt),xx);
    plot(xx,yy,'-','LineWidth',2)
    
end
priorXX = betapdf(xx,priorMix(tt),sum(priorMix)-priorMix(tt));
    plot(xx,priorXX,'--k')
xlim([0 1])
xlabel('Mixing probability')
ylabel('Posterior density')
leg = TypeList; leg{numel(leg)+1} = 'prior';
legend(leg)
hold off
saveas(h,['figures/' modelname '_mixprobs_joint.png'])

h = figure;
hold all
    yy = ksdensity(sum(MP(:,3:4),2),xx);
    plot(xx,yy,'LineWidth',2)
    yy = ksdensity(sum(MP(:,[2 4]),2),xx);
    plot(xx,yy,'LineWidth',2)
    xlim([0 1])
    priorXX = betapdf(xx,sum(priorMix([3 4])),sum(priorMix)-sum(priorMix([3 4])));
    plot(xx,priorXX,'--k')
    legend('A-MP','C-MP','prior')
    xlabel('Mixing probability')
    ylabel('Posterior density')
hold off

saveas(h,['figures/' modelname '_mixprobs_marginal.png'])

%% Mixing export mixing probabilities

eval(['mix_' modelname '= permute(MIXPROBS,[3 2 1]);']);
fpathmix = ['C:\Users\jbland\Documents\LearningSims\' modelname '_mixprobs.mat']
save(fpathmix,['mix_' modelname])


%% Plots - Set some sefault parameters

ci = 0.9; % bayesian credible region to plot

%% individual mixing probabilities -- joint

mMPi = prctile(PRTYPE,0.5,3);
plMPi = prctile(PRTYPE,(1-ci)./2,3);
prMPi = prctile(PRTYPE,1-(1-ci)./2,3);
for tt = 1:4
    h = figure;
    hold all
    data = PRTYPE(:,tt,:);
    plotscript
    
    
    iitype = permute(sum(ITYPE==tt,1),[3 2 1]);
    %F05 = @(x) prctile(cdf('beta',x,iitype,N-iitype),5);
    %F95 = @(x) prctile(cdf('beta',x,iitype,N-iitype),95);
    %l3 = fplot(F05,[0 1],'-b');
    %fplot(F95,[0 1],'-b')
    
    xlabel(['Pr(' TypeList{tt} ')'])
    ylabel('cumulative density')
    
    legend([l1 l2],{'Ind 50th percentile','Ind 5-95th percentile'},'Location','best','Box','off')
    
    
    xlim([-0.05 1.05])
    ylim([0 1])
    whitespace
    hold off
    saveas(h,['figures/' modelname '_mixprobsi_' TypeList{tt} '.png'])
end

%% Just plot some of them
cc = 1;
hsf = figure;
titles = {'(a) Bayes in both tasks','(b) MP in A-task only'};
for tt = [1 3]
    h = subplot(1,2,cc);
    hold all
    title(titles{cc})
    data = PRTYPE(:,tt,:);
    plotscript
    delete(l4)
    
    iitype = permute(sum(ITYPE==tt,1),[3 2 1]);
    %F05 = @(x) prctile(cdf('beta',x,iitype,N-iitype),5);
    %F95 = @(x) prctile(cdf('beta',x,iitype,N-iitype),95);
    %l3 = fplot(F05,[0 1],'-b');
    %fplot(F95,[0 1],'-b')
    
    xlabel(['Pr(' TypeList{tt} ')'])
    ylabel('cumulative density')
    
    if cc==2
        legend([l1 l2],{'Ind 50th percentile','Ind 5th-95th percentile'},'Location','best','Box','off')
    end
    
    xlim([-0.05 1.05])
    ylim([0 1])
    hold off
    cc = cc+1;
end
saveas(hsf,['figures/' modelname '_mixprobsi_2mostCommon.png'])


%% individual mixing probabilities -- marginal

PriMarg = [sum(PRTYPE(:,[3 4],:),2) sum(PRTYPE(:,[2 4],:),2)];
MargList = {'A-MP','C-MP'};
MargType = [(ITYPE == 3 | ITYPE == 4) (ITYPE == 2 | ITYPE == 4)];

for tt = 1:numel(MargList)
    h = figure;
    hold all
    data = PriMarg(:,tt,:);
    plotscript
    
     %F05 = @(x) prctile(cdf('beta',x,sum(MargType(:,tt,:),1),N-sum(MargType(:,tt,:),1)),100*(1-ci)./2);
     %F95 = @(x) prctile(cdf('beta',x,sum(MargType(:,tt,:),1),N-sum(MargType(:,tt,:),1)),100*(1-(1-ci))./2);
     %l3 = fplot(F05,[0 1],'-b');
     %fplot(F95,[0 1],'-b')
    
    xlabel(['Pr(' MargList{tt} ')'])
    ylabel('cumulative density')
    
    legend([l1 l2],'Ind 50th percentile','Ind 5-95th percentile','Location','best')
    
    
    xlim([-0.05 1.05])
    ylim([0 1])
    hold off
    saveas(h,['figures/' modelname '_mixprobsi_marginal_' MargList{tt} '.png'])
    
end

%% Parameters common to all models: gamma and delta

for tt = 1:2
    h = figure;
    hold all
        data  = THETA(:,tt,:);
        plotscript
        
        m = permute(THETA_MU(1,tt,:),[3 1 2]);
        v = permute(THETA_V(tt,tt,:),[3 1 2]);
        F1 = @(x) prctile(normcdf((x-m)./sqrt(v)),100*(1-ci)./2,1);
        F2 = @(x) prctile(normcdf((x-m)./sqrt(v)),100*(1-(1-ci))./2,1);
        
        l3 = fplot(F1,[min(xl) max(xr)],'-b');
        fplot(F2,[min(xl) max(xr)],'-b')
        xlabel(ParamLabels{tt})
        legend([l1 l2 l3],'Ind 50th percentile','Ind 5-95th percentile','Pop cdf: 5th-95th percentile','Location','best')
    hold off
    saveas(h,['figures/' modelname '_params_' ParamTxt{tt} '.png'])
end

%% Parameters common to Bayes model

h = figure;
for pp = 1:numel(PlotThisB)
    p = PlotThisB(pp);
    subplot(numel(PlotThisB)./2,2,pp)
    hold all
        data = THETA(:,p,:);
        m = THETA_MU(1,p,:); v = THETA_V(p,p,:);
        F1 = @(x) prctile(normcdf((x-m(:))./sqrt(v(:))),100*(1-ci)./2);
        F2 = @(x) prctile(normcdf((x-m(:))./sqrt(v(:))),100*(1-(1-ci))./2);
        xx = linspace(min(data(:)),max(data(:)),101)';
        f = zeros(numel(xx),2);
        for kk = 1:numel(xx);
            f(kk,:) = [F1(xx(kk)) F2(xx(kk))];
        end
        eval(['data = ' transform{p} '(data);'])
        eval(['xx = ' transform{p} '(xx);'])
        plotscript
        plot(xx,f,'-b')
        xlim(xlimsB(pp,:))
        xlabel(ParamLabels{p})
    hold off
end
saveas(h,['figures/' modelname '_BayesianParams.png'])

%% Parameters common to MP model


PlotOrder = zeros(numel(PlotThisMP)./2,2); PlotOrder(:) = 1:numel(PlotOrder); PlotOrder = PlotOrder'; PlotOrder= PlotOrder(:);


h = figure;
for pp = 1:numel(PlotThisMP)
    p = PlotThisMP(pp);
    subplot(numel(PlotThisMP)./2,2,pp)
    hold all
        data = THETA(:,p,:);
        m = THETA_MU(1,p,:); v = THETA_V(p,p,:);
        F1 = @(x) prctile(normcdf((x-m(:))./sqrt(v(:))),100*(1-ci)./2);
        F2 = @(x) prctile(normcdf((x-m(:))./sqrt(v(:))),100*(1-(1-ci))./2);
        xx = linspace(min(data(:)),max(data(:)),101)';
        f = zeros(numel(xx),2);
        for kk = 1:numel(xx);
            f(kk,:) = [F1(xx(kk)) F2(xx(kk))];
        end
        eval(['data = ' transform{p} '(data);'])
        eval(['xx = ' transform{p} '(xx);'])
        plotscript
        plot(xx,f,'-b')
        xlim(xlimsMP(pp,:))
        xlabel(['$' ParamLabels{p} '$'])
    hold off
end
saveas(h,['figures/' modelname '_MPParams.png'])









%% MP model - Just alpha

p = PlotThisMP(1);
h = figure;
 hold all
        data = THETA(:,p,:);
        m = THETA_MU(1,p,:); v = THETA_V(p,p,:);
        F1 = @(x) prctile(normcdf((x-m(:))./sqrt(v(:))),100*(1-ci)./2);
        F2 = @(x) prctile(normcdf((x-m(:))./sqrt(v(:))),100*(1-(1-ci))./2);
        xx = linspace(min(data(:)),max(data(:)),101)';
        f = zeros(numel(xx),2);
        for kk = 1:numel(xx);
            f(kk,:) = [F1(xx(kk)) F2(xx(kk))];
        end
        eval(['data = ' transform{p} '(data);'])
        eval(['xx = ' transform{p} '(xx);'])
        plotscript
        plot(xx,f,'-b')
        xlim(xlimsMP(pp,:))
        xlabel(['$' ParamLabels{p} '$'])
        whitespace
    hold off
    saveas(h,['figures/' modelname '_MPalpha.png'])
    

%% New Figure 7
    Apmid = normcdf(THETA(:,7,:));
    Ap0L = Apmid-Apmid.*normcdf(THETA(:,8,:));
    Ap0R = Apmid + (1-Apmid).*normcdf(THETA(:,8,:));
    Apspread = Ap0R-Ap0L;
    Cpmid = normcdf(THETA(:,13,:));
    Cp0L = Cpmid-Cpmid.*normcdf(THETA(:,14,:));
    Cp0R = Cpmid + (1-Cpmid).*normcdf(THETA(:,14,:));
    Cpspread = Cp0R-Cp0L;
Fig7PlotThis = [Apmid Apspread THETA(:,1,:)  exp(THETA(:,2,:))];
xlims        = [0 1   ; 0 1;   -0.5 1.5;         0 40];
xlabels = {'A - midpoint','A - spread','$\gamma$ - risk aversion', '$\lambda$ - choice precision'};
letters = {'(a)','(b)','(c)','(d)'};
hsf = figure('units','inches','Position',[0.1 1 8 8]);

for pp = 1:4
    h = subplot(2,2,pp)
    hold all
    title(letters{pp},'Interpreter','latex')
        data = Fig7PlotThis(:,pp,:);
        plotscript
        FFXX  = linspace(xlims(pp,1),xlims(pp,2),101)';
        F05 = @(x) prctile(mean(Fig7PlotThis(:,pp,:)<=x,1),5);
        F95 = @(x) prctile(mean(Fig7PlotThis(:,pp,:)<=x,1),95);
        FF = zeros(numel(FFXX),2);
        for xx = 1:numel(FFXX)
            FF(xx,:) = [F05(FFXX(xx)) F95(FFXX(xx))];
        end
        
        cdf = plot(FFXX,FF,'-b')
        
        xlim(xlims(pp,:))
        ylim([0 1])
        delete(l4)
        xlabel(xlabels{pp},'Interpreter','latex')
        
        
        if pp==2
            legend([l1 cdf(1)],{'Ind 50th percentile','Pop 5-95th percentile'},'Interpreter','latex','Location','best')
        end
         if pp==4
            legend([l2],{'Ind 5-95th percentile'},'Interpreter','latex','Location','best')
            %whitespace
         end
       
    hold off
end
saveas(hsf,['figures/' modelname '_Fig7.png']) 


    




%% Plot of midpoint of the set of priors in A-task
if strcmp(modelname,'MixBeta')
% nR = 1001;
% pgrid = linspace(0,1,101);
% 
% R = mvnrnd([0 0],[1 0;0 1],nR);
% qL = zeros(SimSize-BurnIn,nR);
% qM = zeros(SimSize-BurnIn,nR);
% qH = zeros(SimSize-BurnIn,nR);
% for ss = 1:(SimSize-BurnIn);
%     q = (repmat(THETA_MU(:,7:8,ss)',[1 nR])+chol(THETA_V(7:8,7:8,ss))*R');
%     qmid = normcdf(q(1,:));
%     qspread   = normcdf(q(2,:));
%     qL(ss,:) = qmid - qmid.*qspread;
%     qH(ss,:) = qmid + (1-qmid).*qspread;
%     qM(ss,:) = qmid;
% end
% QL = zeros(numel(pgrid),2);
% QH = zeros(numel(pgrid),2);
% QM = zeros(numel(pgrid),2);
% for pp = 1:numel(pgrid)
%     p = pgrid(pp);
%     %prctile(data,100*(1-ci)./2,3);
%     QL(pp,1) = prctile(mean(qL<=p,2),100*(1-ci)./2);
%     QL(pp,2) = prctile(mean(qL<=p,2),100-100*(1-ci)./2);
%     QH(pp,1) = prctile(mean(qH<=p,2),100*(1-ci)./2);
%     QH(pp,2) = prctile(mean(qH<=p,2),100-100*(1-ci)./2);
%     QM(pp,2) = prctile(mean(qM<=p,2),100-100*(1-ci)./2);
%     QM(pp,2) = prctile(mean(qM<=p,2),100-100*(1-ci)./2);
% end




    pmid = normcdf(THETA(:,7,:));
    p0L = pmid-pmid.*normcdf(THETA(:,8,:));
    p0R = pmid + (1-pmid).*normcdf(THETA(:,8,:));
    data = 0.5.*(p0L+p0R);
    ci = 0.95;
    h = figure;

    hold all
        plotscript
        plot([0.5 0.5],[0 1],'--k')
        xlim([0 1])
        ylim([0 1])
        CP = prctile(data,100*[0.025 0.975],3);
        midIn = CP(:,1,1)<0.5 & CP(:,1,2)>0.5;
        disp('Fraction of subjects whose 95% cr edible region')
        disp('of the midpoint of their set of priors covers p=0.5:')
        disp([num2str(mean(midIn))])
        
        F05 = @(x) normcdf(prctile((norminv(x)-THETA_MU(1,7,:))./sqrt(THETA_V(7,7,:)),100-100*(1-ci)./2));
        F95 = @(x) normcdf(prctile((norminv(x)-THETA_MU(1,7,:))./sqrt(THETA_V(7,7,:)),100*(1-ci)./2));
        fplot(F05,[0 1],'-b');
        l3 = fplot(F95,[0 1],'-b');
        
        %l3 = plot(pgrid,QM,'-b');
        xlabel('Midpoint in set of priors')
        legend([l1 l2 l3],'Ind 50th percentile','Ind 5-95th percentile','Pop 5-95th percentile','Location','best')
    hold off
    saveas(h,['figures/' modelname '_MP_midpoint.png'])


    
end




 



%% Heatmaps of related parameters



gs = 101;
for kk = 1:size(PlotThis1,1);
    m = THETA_MU(1,PlotThis1(kk,:),:);
    V = THETA_V(PlotThis1(kk,:),PlotThis1(kk,:),:);
    xx = linspace(PlotRange1(kk,1),PlotRange1(kk,2),gs);
    [AA,CC] = meshgrid(xx,xx);
    F = zeros(size(AA));
    
        f = zeros([numel(AA),size(m,3)]);
        for ss = 1:size(m,3)
           f(:,ss) = exp(logmvnpdf([AA(:) CC(:)],m(:,:,ss),V(:,:,ss))) ;
        end
        mf = median(f,2);
        F(:) = mf(:);
       
    eval(['AA = ' transform{PlotThis1(kk,1)} '(AA);'])
    eval(['CC = ' transform{PlotThis1(kk,2)} '(CC);'])
        
    h = figure;
    mesh(AA,CC,F)
    hold all
        xlabel(ParamLabels{PlotThis1(kk,1)})
        ylabel(ParamLabels{PlotThis1(kk,2)})
    hold off
    
    saveas(h,['figures/' modelname '_2way' ParamTxt{PlotThis1(kk,1)} '.png'])
end

%% "regression" table

test = 0.05;


fid = fopen(['figures/' modelname '.txt'],'w');

fprintf(fid,'\\begin{tabular}{l %s}\\hline\\hline \n',repmat('r',[1 numel(ParamTable)]));

fprintf(fid,'&$%s$',ParamTable{:});
fprintf(fid,'\\\\ \n transform');
fprintf(fid,'& %s',transform{:});
fprintf(fid,'\\\\ \\hline \n ');

fprintf(fid,'{\\sc Mean}');
xm = mean(THETA_MU,3);
xs = std(THETA_MU,[],3);
fprintf(fid,'&%4.4f',xm);
fprintf(fid,'\\\\ \n ');
for kk = 1:numel(ParamTable)
    ci = prctile(data,100.*[test/2 1-test/2],3);
            if 0<ci(1) || 0>ci(2)
                star = '*';
            else
                star = '';
            end
    fprintf(fid,'&(%2.2f)%s',xs(kk),star);
end
fprintf(fid,'\\\\ \\hline \n ');

fprintf(fid,'{\\sc Variance} \\\\ \n');
for kk = 1:numel(ParamTable)
    fprintf(fid,'$%s$',ParamTable{kk});
    for cc = 1:numel(ParamTable)
        if kk>=cc
            fprintf(fid,'&%4.4f',mean(THETA_V(kk,cc,:),3));
        else
            fprintf(fid,'&\\multicolumn{1}{c}{-}');
        end
    end
    fprintf(fid,'\\\\ \n ');
    for cc = 1:numel(ParamTable)
        if kk>=cc
            ci = prctile(THETA_V(kk,cc,:),100.*[test/2 1-test/2],3);
            if (0<ci(1) || 0>ci(2)) && kk~=cc
                star = '*';
            elseif kk==cc
                star = '$^a$';
            else
                star = '';
            end
            fprintf(fid,'&(%2.2f)%s',std(THETA_V(kk,cc,:),[],3),star);
        else
            %fprintf(fid,'&\\multicolumn{1}{c}{-}');
        end
    end
    fprintf(fid,'\\\\ \n ');
end


fprintf(fid,' \\hline\\hline \n ');
fprintf(fid,'\\multicolumn{%s}{l}{Table shows posterior means (standard deviations)} \\\\ \n', int2str(numel(ParamTable)+1));
fprintf(fid,'\\multicolumn{%s}{l}{* indicates that a %s percent Bayesian credible region does not include zero} \\\\ \n', int2str(numel(ParamTable)+1),int2str(100*(1-test)));
fprintf(fid,'\\multicolumn{%s}{l}{$^a$ Stars are supressed because these parameters can only be positive} \n', int2str(numel(ParamTable)+1));

fprintf(fid,'\\end{tabular} \n');
fclose(fid);

%% Actual parameter table - population distribution
% i.e. the above table gives the 
actualPop = 1; % set ~=1 to skip this
if actualPop ==1
    rs = 1000;
    SimParamsPop = zeros(rs,14,SimSize-BurnIn);
    randn('state',42);
    R = randn(rs,14);
    ActualLabels = {'$\gamma$','$\log\lambda$','$\log N_0$','$p_0$','$\alpha$','$\log N_0$','$p_{0mid}$','$p_{0sp}$','$\log N_0$','$p_0$','$\alpha$','$\log N_0$','$p_{0mid}$','$p_{0sp}$'};
    ActualSIM = zeros(rs,14,SimSize-BurnIn);
    ActualSIM_mean = zeros(1,14,SimSize-BurnIn);
    ActualSIM_median = zeros(1,14,SimSize-BurnIn);
    ActualSIM_cov = zeros(14,14,SimSize-BurnIn);
    ActualSIM_corr = zeros(14,14,SimSize-BurnIn);
    for ss = 1:(SimSize-BurnIn)
        % This gets rid of some rounding errors in the 
        TV = (THETA_V(:,:,ss)+THETA_V(:,:,ss)')./2;
        
        SimParamsPop = (repmat(THETA_MU(:,:,ss)',[1 rs])+chol(TV)*R')';
        ActualParams = SimParamsPop;
        % gamma: do nothing
        % N_0s: e^X
        %ActualParams(:,[3 6 9 12]) = exp(SimParamsPop(:,[3 6 9 12]));
        % priors for Bayesians, alpha for MPs
        ActualParams(:,[4 5 10 11]) = normcdf(SimParamsPop(:,[4 5 10 11]));
        % set of priors
        % A-task
        px = normcdf(SimParamsPop(:,7,:));
        p0L = px-px.*normcdf(SimParamsPop(:,8,:));
        p0R = px + (1-px).*normcdf(SimParamsPop(:,8,:));
        ActualParams(:,7) = 0.5.*(p0L+p0R);
        ActualParams(:,8) = (p0R-p0L);
        % C-task
        px = normcdf(SimParamsPop(:,13,:));
        p0L = px-px.*normcdf(SimParamsPop(:,14,:));
        p0R = px + (1-px).*normcdf(SimParamsPop(:,14,:));
        ActualParams(:,13) = 0.5.*(p0L+p0R);
        ActualParams(:,14) = (p0R-p0L);
        ActualSIM(:,:,ss) = ActualParams;
        
        ActualSIM_mean(:,:,ss) = mean(ActualSIM(:,:,ss),1);
        ActualSIM_median(:,:,ss) =  median(ActualSIM(:,:,ss),1);
        ActualSIM_cov(:,:,ss) = cov(ActualSIM(:,:,ss),1);
        ActualSIM_corr(:,:,ss) = corr(ActualSIM(:,:,ss));
    
        disp([ss./(SimSize-BurnIn)])
    end
    
end

%%

VCORR = ActualSIM_corr;
scaling = [0 0 0  0 0 0 0 0 0 0 0 0 0 0];
for vv = 1:14
    VCORR(vv,vv,:) = sqrt(ActualSIM_cov(vv,vv,:));
    if scaling(vv)~=0
        VCORR(vv,vv,:) = exp(log(sqrt(ActualSIM_cov(vv,vv,:)))-scaling(vv)*log(10));
    end
end



test = 0.05;



%% all parameters in one table
fid = fopen(['figures/' modelname '_TRANSFORMED.txt'],'w');

fprintf(fid,'\\begin{tabular}{l %s}\\hline\\hline \n',repmat('r',[1 numel(ParamTable)]));
fprintf(fid,'&\\multicolumn{2}{c}{{\\sc Common}}');
fprintf(fid,'&\\multicolumn{2}{c}{{\\sc Ambiguous - Bayes}}');
fprintf(fid,'&\\multicolumn{4}{c}{{\\sc Ambiguous - Multiple Priors}}');
fprintf(fid,'&\\multicolumn{2}{c}{{\\sc Compound - Bayes}}');
fprintf(fid,'&\\multicolumn{4}{c}{{\\sc Compound - Multiple Priors}}');
fprintf(fid,'\\\\ \n')
for kk = 1:numel(ActualLabels)
    fprintf(fid,'&\\multicolumn{1}{c}{%s',ActualLabels{kk});
    if scaling(kk)~=0
        fprintf(fid,' $\\times 10^{-%2.0f}$',scaling(kk));
    end
    fprintf(fid,'}');
end
%fprintf(fid,'\\\\ \n transform');
%fprintf(fid,'& %s',transform{:});
fprintf(fid,'\\\\ \\hline \n ');

  fprintf(fid,'\\multicolumn{15}{l}{{\\sc Mean}} \\\\');
  data = ActualSIM_mean;
   xm = mean(data,3)./10.^scaling;
   xs = std(data,[],3)./10.^scaling;
   fprintf(fid,'&%1.2f',xm);
   fprintf(fid,'\\\\ \n ');
   for kk = 1:numel(ParamTable)
       ci = prctile(data,100.*[test/2 1-test/2],3);
               if 0<ci(1) || 0>ci(2)
                   star = '*';
               else
                   star = '';
               end
               if kk~=1
                 star = '$^a$';
               end
       fprintf(fid,'&(%1.2f)%s',xs(kk),star);
       %fprintf(fid,'&[%2.2f, %2.2f]',ci(1),ci(2));
   end
  fprintf(fid,'\\\\ \n ');

fprintf(fid,'\\multicolumn{15}{l}{{\\sc Median} ~ (unscaled)} \\\\');
data = ActualSIM_median;
xm = mean(data,3);
xs = std(data,[],3);
fprintf(fid,'&%2.2f',xm);
fprintf(fid,'\\\\ \n ');
for kk = 1:numel(ParamTable)
    ci = prctile(data,100.*[test/2 1-test/2],3);
            if 0<ci(1) || 0>ci(2)
                star = '*';
            else
                star = '';
            end
            if kk~=1
                star = '$^a$';
            end
    fprintf(fid,'&(%2.2f)%s',xs(kk),star);
    %fprintf(fid,'&[%2.2f, %2.2f]',ci(1),ci(2));
end
fprintf(fid,'\\\\ \\hline \n ');

 fprintf(fid,'\\multicolumn{15}{l}{{\\sc Variance \\& Correlation}} \\\\ \n');
  for kk = 1:numel(ParamTable)
      fprintf(fid,'%s',ActualLabels{kk});
     for cc = 1:numel(ParamTable)
         if kk>=cc
             fprintf(fid,'&%1.2f',mean(VCORR(kk,cc,:),3));
             %fprintf(fid,'&%1.2f',median(VCORR(kk,cc,:),3));
         else
             fprintf(fid,'&\\multicolumn{1}{c}{-}');
         end
     end
     fprintf(fid,'\\\\ \n ');
     for cc = 1:numel(ParamTable)
         if kk>=cc
             ci = prctile(VCORR(kk,cc,:),100.*[test/2 1-test/2],3);
             if (0<ci(1) || 0>ci(2)) && kk~=cc
                 star = '*';
             elseif kk==cc
                 star = '$^a$';
             else
                 star = '';
             end
             if kk==cc
                 fprintf(fid,'&[%2.2f, %2.2f]',ci(1),ci(2));
             else
                 fprintf(fid,'&(%2.2f)%s',std(VCORR(kk,cc,:),[],3),star);
             end
             
         else
             %fprintf(fid,'&\\multicolumn{1}{c}{-}');
         end
     end
     fprintf(fid,'\\\\ \n ');
 end


fprintf(fid,' \\hline\\hline \n ');
fprintf(fid,'\\multicolumn{%s}{l}{Table shows posterior means (standard deviations)} \\\\ \n', int2str(numel(ParamTable)+1));
fprintf(fid,'\\multicolumn{%s}{l}{* indicates that a %s percent Bayesian credible region does not include zero} \\\\ \n', int2str(numel(ParamTable)+1),int2str(100*(1-test)));
fprintf(fid,'\\multicolumn{%s}{l}{$^a$ Stars are supressed because these parameters can only be positive} \n', int2str(numel(ParamTable)+1));

fprintf(fid,'\\end{tabular} \n');
fclose(fid);

%% all parameters in one table, but cutting out the CIs

ShowStars = [1 0 1 0 0 1 0 0 1 0 0 1 0 0];

ordering = [1 2     ... % common
            9 10    ... % compound BayesMixProb
            3 4     ... % ambiguous Bayes
            5 6 7 8];   % ambiguous MP

fid = fopen(['figures/' modelname '_TRANSFORMED_smaller.txt'],'w');

fprintf(fid,'\\begin{tabular}{l %s}\\hline\\hline \n',repmat('r',[1 numel(ordering)]));
fprintf(fid,'&\\multicolumn{2}{c}{{\\sc Common}}');
fprintf(fid,'&\\multicolumn{2}{c}{{\\sc Compound - Bayes}}');
fprintf(fid,'&\\multicolumn{2}{c}{{\\sc Ambiguous - Bayes}}');
fprintf(fid,'&\\multicolumn{4}{c}{{\\sc Ambiguous - Multiple Priors}}');

%fprintf(fid,'&\\multicolumn{4}{c}{{\\sc Compound - Multiple Priors}}');
fprintf(fid,'\\\\ \n')
for jj = 1:numel(ordering)
    kk = ordering(jj);
    fprintf(fid,'&\\multicolumn{1}{c}{%s',ActualLabels{kk});
    if scaling(kk)~=0
        fprintf(fid,' $\\times 10^{-%2.0f}$',scaling(kk));
    end
    fprintf(fid,'}');
end
%fprintf(fid,'\\\\ \n transform');
%fprintf(fid,'& %s',transform{:});
fprintf(fid,'\\\\ \\hline \n ');

  fprintf(fid,'\\multicolumn{ %s}{l}{{\\sc Mean}} \\\\',int2str(numel(ordering)));
  data = ActualSIM_mean;
   xm = mean(data,3)./10.^scaling;
   xs = std(data,[],3)./10.^scaling;
   fprintf(fid,'&%1.2f',xm(ordering));
   fprintf(fid,'\\\\ \n ');
   for jj = 1:numel(ordering)
       kk = ordering(jj);
       ci = prctile(data,100.*[test/2 1-test/2],3);
               if 0<ci(1,kk,1) || 0>ci(1,kk,2)
                   star = '*';
               else
                   star = '';
               end
               if ShowStars(kk)==0
                 star = '$^a$';
               end
       fprintf(fid,'&(%1.2f)%s',xs(kk),star);
       %fprintf(fid,'&[%2.2f, %2.2f]',ci(1),ci(2));
   end
   
   fprintf(fid,'\\\\ \\hline \n ');
 
 fprintf(fid,'\\multicolumn{5}{l}{{\\sc Variance \\& Correlation}} \\\\ \n');
  for ii = 1:numel(ordering)
      kk = ordering(ii);
      fprintf(fid,'%s',ActualLabels{kk});
     for jj = 1:numel(ordering)
         cc = ordering(jj);
         ci = prctile(VCORR(kk,cc,:),100.*[test/2 1-test/2],3);
             if (0<ci(1) || 0>ci(2)) && kk~=cc
                 star = '*';
             elseif kk==cc
                 star = '$^a$';
             else
                 star = '';
             end
         if ii>=jj
             fprintf(fid,'&%1.2f%s',mean(VCORR(kk,cc,:),3),star);
             %fprintf(fid,'&%1.2f',median(VCORR(kk,cc,:),3));
         else
             fprintf(fid,'&\\multicolumn{1}{c}{~}');
         end
     end
     fprintf(fid,'\\\\ \n ');
 end


fprintf(fid,' \\hline\\hline \n ');
fprintf(fid,'\\multicolumn{%s}{l}{Table shows posterior means (standard deviations)} \\\\ \n', int2str(numel(ordering)+1));
fprintf(fid,'\\multicolumn{%s}{l}{* indicates that a %s percent Bayesian credible region does not include zero} \\\\ \n', int2str(numel(ordering)+1),int2str(100*(1-test)));
fprintf(fid,'\\multicolumn{%s}{l}{$^a$ Stars are supressed because these parameters can only be positive} \n', int2str(numel(ordering)+1));

fprintf(fid,'\\end{tabular} \n');
fclose(fid);


%% Correlation/means table - smaller

TabThese = [1 2 5:8];

%D = ActualSIM(:,[1 2 5 7],:);
D = ActualSIM(:,TabThese,:);
test = 0.05;
mD = mean(D,1);

 labels = ActualLabels(TabThese);

fid = fopen(['figures/' modelname '_SmallCorrelations.txt'],'w');

fprintf(fid,'\\begin{tabular}{l %s}\\hline\\hline \n',repmat('r',[1 numel(labels)]));

fprintf(fid,'&\\multicolumn{1}{c}{%s}',labels{:});
fprintf(fid,'\\\\ \\hline \n ');

fprintf(fid,'{\\sc Mean}');
xm = mean(mD,3);
xs = std(mD,[],3);
fprintf(fid,'&%1.2f',xm); %fprintf(fid,'&%4.4f',xm);
fprintf(fid,'\\\\ \n ');
for kk = 1:numel(labels)
    ci = prctile(data,100.*[test/2 1-test/2],3);
            if 0<ci(1) || 0>ci(2)
                star = '*';
            else
                star = '';
            end
    fprintf(fid,'&(%1.2f)%s',xs(kk),star);
end
fprintf(fid,'\\\\ \\hline \n ');

%cD = zeros(numel(labels),numel(labels),size(D,3));
%vD = zeros(numel(labels),numel(labels),size(D,3));
xD = zeros(numel(labels),numel(labels),size(D,3));
for ss = 1:size(D,3)
    cD = corr(D(:,:,ss));
    vD = var(D(:,:,ss));
    
    xD(:,:,ss) = cD-eye(numel(labels))+diag(vD);
end



fprintf(fid,'\\multicolumn{3}{l}{{\\sc Variance / correlation}} \\\\ \n');
for kk = 1:numel(labels)
    fprintf(fid,'%s',labels{kk});
     for cc = 1:numel(labels)
        if kk>=cc
            fprintf(fid,'&%1.2f',mean(xD(kk,cc,:),3));
        else
            fprintf(fid,'&\\multicolumn{1}{c}{-}');
        end
    end
    
    fprintf(fid,'\\\\ \n ');
    for cc = 1:numel(labels)
        if kk>=cc
            ci = prctile(xD(kk,cc,:),100.*[test/2 1-test/2],3);
            if (0<ci(1) || 0>ci(2)) && kk~=cc
                star = '*';
            elseif kk==cc
                star = '$^a$';
            else
                star = '';
            end
            fprintf(fid,'&(%1.2f)%s',std(xD(kk,cc,:),[],3),star);
        else
            %fprintf(fid,'&\\multicolumn{1}{c}{-}');
        end
    end
    fprintf(fid,'\\\\ \n ');
end


fprintf(fid,' \\hline\\hline \n ');
fprintf(fid,'\\multicolumn{%s}{p{0.7\\textwidth}}{Table shows posterior means (standard deviations)} \\\\ \n', int2str(numel(labels)+1));
fprintf(fid,'\\multicolumn{%s}{p{0.7\\textwidth}}{* indicates that a %s percent Bayesian credible region does not include zero} \\\\ \n', int2str(numel(labels)+1),int2str(100*(1-test)));
fprintf(fid,'\\multicolumn{%s}{p{0.7\\textwidth}}{$^a$ Stars are supressed because these parameters can only be positive} \n', int2str(numel(labels)+1));

fprintf(fid,'\\end{tabular} \n');
fclose(fid);


%% Heatmap of p_mid and N_0 from A-task
set(0,'defaulttextInterpreter','latex') 
plot_p = ActualSIM(:,7,:); plot_p = plot_p(:);
plot_N = ActualSIM(:,6,:); plot_N = plot_N(:);

keep = plot_N<=400;

plot_p = plot_p(keep);
plot_N = plot_N(keep);

PP = repmat(0:0.01:1,[101 1]);
NN = repmat(0:0.1:10,[101 1])';

plotpoints = median([plot_p plot_N]);
h = figure;
hold all
    
    %ksdensity([plot_p plot_N],[PP(:) NN(:)]);
    ksdensity([plot_p plot_N])
    colormap(jet)
    for kk = 1:size(plotpoints,1)
        %annotation('textarrow',[plotpoints(kk,1), plotpoints(kk,1)+0.1],[plotpoints(kk,2) plotpoints(kk,2)+0.1],'String','y = x ')
        %text(,'A','HorizontalAlignment','right')
        %annotation('textbox',plotpoints(kk,1),plotpoints(kk,2))
        plot(plotpoints(kk,1),plotpoints(kk,2),'or')
    end
    %heatmap(XX,YY);
    xlabel('$p_{mid}$')
    %ylim([0 1000])
    ylabel('$\log N_0$')
    xlim([0 1])
    ylim([0 400])
hold off
saveas(h,['figures/' modelname '_heatmap_N0p0mid.png'])
    
%% Heatmap of p_mid and alpha from A-task
set(0,'defaulttextInterpreter','latex') 
plot_p = ActualSIM(:,7,:); plot_p = plot_p(:);
plot_N = ActualSIM(:,5,:); plot_N = plot_N(:);

%keep = plot_N<=400;

%plot_p = plot_p(keep);
%plot_N = plot_N(keep);

PP = repmat(0:0.01:1,[101 1]);
NN = repmat(0:0.1:10,[101 1])';

plotpoints = median([plot_p plot_N]);
h = figure;
hold all
    
    %ksdensity([plot_p plot_N],[PP(:) NN(:)]);
    ksdensity([plot_p plot_N])
    colormap(jet)
    for kk = 1:size(plotpoints,1)
        %annotation('textarrow',[plotpoints(kk,1), plotpoints(kk,1)+0.1],[plotpoints(kk,2) plotpoints(kk,2)+0.1],'String','y = x ')
        %text(,'A','HorizontalAlignment','right')
        %annotation('textbox',plotpoints(kk,1),plotpoints(kk,2))
        plot(plotpoints(kk,1),plotpoints(kk,2),'or')
    end
    %heatmap(XX,YY);
    xlabel('$p_{mid}$')
    %ylim([0 1000])
    ylabel('$\alpha$')
    xlim([0 1])
    ylim([0 1])
hold off
saveas(h,['figures/' modelname '_heatmap_alphap0mid.png'])
   


%% fraction not bayes

NBfunA = @(x) NotBayes_MPbeta([E1_Choice_A E1_Choice_A2],[E1_Prob_A E1_Prob_A2],[E1_SafeOption_A E1_SafeOption_A2],[x(:,1:2) x(:,5:8)],[E1_NumDraws_A E1_NumDraws_A2],[E1_Draw_A E1_Draw_A2],[E1_DecisionColor_A E1_DecisionColor_A2]);
NBfunC = @(x) NotBayes_MPbeta(E1_Choice_C,E1_Prob_C,E1_SafeOption_C,[x(:,1:2) x(:,11:14)],E1_NumDraws_C,E1_Draw_C,E1_DecisionColor_C);
NOTBAYES_A = zeros(N,size(THETA,3));
NOTBAYES_C = zeros(N,size(THETA,3));
for ss = 1:size(THETA,3)
    NOTBAYES_A(:,ss) = sum(PRTYPE(:,3:4),2).*NBfunA(THETA(:,:,ss));
     NOTBAYES_C(:,ss) = sum(PRTYPE(:,[2 4]),2).*NBfunC(THETA(:,:,ss));
    disp(ss/size(THETA,3))
end

%% plot fractions not bayes

notbayesA = [mean(NOTBAYES_A,2) prctile(NOTBAYES_A,[5 95],2)];
notbayesC = [mean(NOTBAYES_C,2) prctile(NOTBAYES_C,[5 95],2)];

ci = 0.05;
h = figure;
h1 = subplot(1,2,1)
hold all
    title('(a) A-task')
    data = permute(NOTBAYES_A,[1 3 2]);
    xm = prctile(data,50,3);
    xl = prctile(data,100*(1-ci)./2,3);
    xr = prctile(data,100-100*(1-ci)./2,3);
    [xs,I] = sort(xm);
    %for ii = 1:numel(xm)
        l1 = plot(xm(I),(1:numel(xm))'./numel(xm),'.k');
    for ii = 1:numel(xm)
        l2 = plot([xl(I(ii)),xr(I(ii))],repmat(ii./numel(xm),[1 2]),'-k');
    end
    ylabel('cumulative density')
    %xlabel('Fraction of potentially different decisions')
    %legend([l1 l2],{'Ind 50th percentile','Ind 5-95th percentile'},'Location','best', 'Box','off')
    xlim([0 1])
    ks = [];
    
    %saveas(h,['figures/' modelname '_NotBayesA.png'])

hold off

h2 = subplot(1,2,2)
hold all
    title('(b) C-task')
    data = permute(NOTBAYES_C,[1 3 2]);
    xm = prctile(data,50,3);
    xl = prctile(data,100*(1-ci)./2,3);
    xr = prctile(data,100-100*(1-ci)./2,3);
    [xs,I] = sort(xm);
    %for ii = 1:numel(xm)
        l1 = plot(xm(I),(1:numel(xm))'./numel(xm),'.k');
    for ii = 1:numel(xm)
        l2 = plot([xl(I(ii)),xr(I(ii))],repmat(ii./numel(xm),[1 2]),'-k');
    end
    %ylabel('cumulative density')
    %xlabel('Fraction of potentially different decisions')
    legend([l1 l2],{'Ind 50th percentile','Ind 5th-95th percentile'},'Location','best', 'Box','off','Interpreter','latex')
    xlim([0 1])
hold off
p1=get(h1,'position');
p2=get(h2,'position');
p3=get(h1,'position');
p4=get(h2,'position');
height=p1(2)+p1(4)-p4(2);
width=p4(1)+p4(3)-p3(1);
h5=axes('position',[p3(1) p3(2) width height],'visible','off'); 
h5.XLabel.Visible='on'
h5.YLabel.Visible='on'
axes(h5)
%ylabel('test')
xlabel('Fraction of potentially different decisions')
saveas(h,['figures/' modelname '_NotBayes.png'])

% for tt = 1:2
%     h = figure;
%     hold all
%         data  = THETA(:,tt,:);
%         plotscript
%         
%         m = permute(THETA_MU(1,tt,:),[3 1 2]);
%         v = permute(THETA_V(tt,tt,:),[3 1 2]);
%         F1 = @(x) prctile(normcdf((x-m)./sqrt(v)),100*(1-ci)./2,1);
%         F2 = @(x) prctile(normcdf((x-m)./sqrt(v)),100*(1-(1-ci))./2,1);
%         
%         l3 = fplot(F1,[min(xl) max(xr)],'-b');
%         fplot(F2,[min(xl) max(xr)],'-b')
%         xlabel(ParamLabels{tt})
%         legend([l1 l2 l3],'Ind 50th percentile','Ind 5-95th percentile','Pop cdf: 5th-95th percentile','Location','best')
%     hold off
%     saveas(h,['figures/' modelname '_params_' ParamTxt{tt} '.png'])
% end

%% Collect some things needed to compare models



LIKE = zeros(N,4,size(MIXPROBS,3));
LIKEparams = zeros(N,1,size(MIXPROBS,3));
for kk = 1:4
    eval(['ll_' int2str(kk) '=ll{' int2str(kk) '};'])
end

for ss = 1:size(MIXPROBS,3)
    for kk = 1:4
        eval(['LIKE(:,kk,ss) = ll_' int2str(kk) '(THETA(:,:,ss));'])
        LIKEparams(:,1,ss) = logmvnpdf(THETA(:,:,ss),THETA_MU(:,:,ss),THETA_V(:,:,ss))';
    end
    

    disp(ss./size(MIXPROBS,3))
end



 
 MARGLIKE = LIKEparams;
 
 for ii = 1:size(MARGLIKE,1)
        litmp = permute(LIKE(ii,:,:),[3,2,1]);
        itypetmp =permute(ITYPE(ii,:,:),[3,2,1]);
        liketmp = zeros(1,size(MARGLIKE,3));
        for tt=1:4
            liketmp(itypetmp==tt) = litmp(itypetmp==tt,tt);
        end
        
        MARGLIKE(ii,1,:) = MARGLIKE(ii,1,:) + permute(liketmp,[3,1,2]);
 end
     
 
 keeplist = 'MARGLIKE MIXPROBS THETA';
 
 eval(['save C:\Users\jbland\Documents\LearningSims\modelcomp_' modelname '.mat ' keeplist ])
 
%%

close all
