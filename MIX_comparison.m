%% Comparison of mixture models
clear; clc

TypeList = {'A-B C-B','A-B C-MP','A-MP C-B','A-MP C-MP'};

mixlist = {'MixSimplexGeneral','MixBeta','MixBetaRestricted','MixBeta_Both'};
keeplist = {'MARGLIKE','MIXPROBS','THETA'};

for mm = 1:numel(mixlist)
    eval(['load C:\Users\jbland\Documents\LearningSims\modelcomp_' mixlist{mm} '.mat'])
    
    for kk = 1:numel(keeplist)
        eval([keeplist{kk} '_' int2str(mm) '=' keeplist{kk} ';'])
        eval(['clear ' keeplist{kk}  ])
        if strcmp(keeplist{kk},"MARGLIKE")
            eval(['ml = MARGLIKE_' int2str(mm) ';'])
            
            sprintf([mixlist{mm} ' marginal likelihood  = ' num2str(mean(sum(ml,1),3))])
        end
    end
    
end
kill

%%
ML1 = permute(sum(MARGLIKE_1,1),[3 2 1]);
ML2 = permute(sum(MARGLIKE_2,1),[3 2 1]);
ML3 = permute(sum(MARGLIKE_3,1),[3 2 1]);
ML4 = permute(sum(MARGLIKE_4,1),[3 2 1]);

minsize = min([numel(ML1),numel(ML2),numel(ML3),numel(ML4)])

%ksdensity(ML2-ML1)
DML = ML2-ML1;
h = figure;
hold all
plot(DML,'.k')
title(['Beta - Simplex:' 'mean log posterior odds ratio = ' num2str(mean(DML))])
xlabel('Simulation step (after burn-in)')
ylabel('Log posterior odds ratio')
hold off
saveas(h,'figures/MIXCOMPARE.png')

%ksdensity(ML2-ML1)
DML = ML2-ML4;
h = figure;
hold all
plot(DML,'.k')
title(['Beta - Restricted Beta:' 'mean log posterior odds ratio = ' num2str(mean(DML))])
xlabel('Simulation step (after burn-in)')
ylabel('Log posterior odds ratio')
hold off
saveas(h,'figures/MIXCOMPARE_BetaRestriction.png')

%% Bayes Factor

    addThis = -min([ML1;ML2]);
    MLexp1 = mean(exp(ML1+addThis));
    MLexp2 = mean(exp(ML2+addThis));
    log10BayesFactorBS = log10(MLexp2)-log10(MLexp1);
    
    addThis = -max([ML2;ML4]);
    MLexp1 = mean(exp(ML4+addThis));
    MLexp2 = mean(exp(ML2+addThis));
    log10BayesFactorBBr = log10(MLexp2)-log10(MLexp1);
    
    ["Beta - Simplex" num2str(log10BayesFactorBS);
      "Beta - Beta restricted" num2str(log10BayesFactorBBr)  ]
    
%%

BayesFactor = MLexp1./MLexp2

1./BayesFactor


%%
p2 = 1./(1+exp(-DML));

MIXPROBS_W = permute(MIXPROBS_2,[3 2 1]).*repmat(p2,[1 4])+permute(MIXPROBS_1,[3 2 1]).*repmat(1-p2,[1 4]);

%% Problem: big shift in marginal likelihood of MixBeta

h = figure;
for tt = 1:4
   subplot(2,2,tt)
   hold all
        title(TypeList(tt))
        ksdensity(permute(MIXPROBS_1(1,tt,:),[3 1 2]))
        ksdensity(permute(MIXPROBS_2(1,tt,:),[3 1 2]))
        [F,X] = ksdensity(MIXPROBS_W(:,tt));
        plot(X,F,'--k')
        xlim([0 1])
        if tt==4
            legend('Simplex','Beta','BMA')
        end
   hold off
end
saveas(h,'figures/MIXCOMPARE_BMA.png')