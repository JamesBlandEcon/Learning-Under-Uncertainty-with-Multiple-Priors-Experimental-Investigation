%% simulates the posterior distribution assuming that:
% either MP of subjective bayes 
% parameters and type can be different between tasks
% Beta priors used
clear; clc;

load Exp1.mat
load Expt2.mat

Experiment = 1;

doML = 0; % set ==1 to do individual ML estimates (slow and ultimately useless when there are a lot of parameters)

E1_DecisionColor_A= strcmp(E1_DecisionColor_A,'b');
E1_DecisionColor_R= strcmp(E1_DecisionColor_R,'b');
E1_DecisionColor_A2= strcmp(E1_DecisionColor_A2,'b');
E1_DecisionColor_C= strcmp(E1_DecisionColor_C,'b');

CombineList = {"Choice","Prob","SafeOption","NumDraws","Draw","DecisionColor"};

for cc = 1:length(CombineList)
    
    c = CombineList{cc};

    str = strcat('X= [E1_',c,'_A E1_',c,'_A2];');
    eval(str)
    eval(strcat('AddThis=A_',c,';'))
    Y = [X; AddThis NaN*zeros(size(AddThis,1),size(X,2)-size(AddThis,2))];
    eval(strcat('A',c,'=Y;'))
    
    for tt = ["C" "R"]
        str = strcat('X=E1_',c,"_",tt,';');
        eval(str)
        eval(strcat('AddThis=',tt,'_',c,';'))
        Y = [X; AddThis NaN*zeros(size(AddThis,1),size(X,2)-size(AddThis,2))];
        eval(strcat(tt,c,'=Y;'))
    end

end






doML = 0; % set ==1 to do individual ML estimates (slow and ultimately useless when there are a lot of parameters)




modelname = 'MixBeta_Both'



%% Set up the likelihood function


%parameter vector theta is:
%    parameter   task     type starting_value
% ------------------------------------------------
%01  gamma       all      all       0.5
%02  lambda      all      all       exp(2)
%03  N0          A        Bayes     2  % i.e. uniform prior
%04  p0          A        Bayes     0.5
%05  alpha       A        MP
%06  N0          A        MP
%07  p0mid       A        MP
%08  p0spread    A        MP
%09  N0          C        Bayes     2
%10  p0          C        Bayes     0.5
%11  alpha       C        MP
%12  N0          C        MP
%13  p0mid       C        MP
%14  p0spread    C        MP

[N,T] = size(AChoice);

theta0 = repmat([0.5 2 log(2) 0 0 0 0 0 log(2) 0 0 0 0 0],[N 1]);


% 00 - Bayes in both tasks
ll{1} = @(x)    A_loglike_BAYES_Beta(AChoice,ASafeOption,x(:,1:4),ANumDraws,ADraw,ADecisionColor)...
              +A_loglike_BAYES_Beta(CChoice,CSafeOption,[x(:,1:2) x(:,9:10)],CNumDraws,CDraw,CDecisionColor)...
              +R_loglike2(RChoice,RProb,RSafeOption,x(:,1:2));
% 01 - A-Bayes & C-MP
ll{2} = @(x)    A_loglike_BAYES_Beta(AChoice,ASafeOption,x(:,1:4),ANumDraws,ADraw,ADecisionColor)...
              +A_loglike_MPbeta_black(CChoice,CProb,CSafeOption,[x(:,1:2) x(:,11:14)],CNumDraws,CDraw,CDecisionColor)...
              +R_loglike2(RChoice,RProb,RSafeOption,x(:,1:2));

% 10 - A-MP & C-Bayes
ll{3} = @(x)    A_loglike_MPbeta_black(AChoice,AProb,ASafeOption,[x(:,1:2) x(:,5:8)],ANumDraws,ADraw,ADecisionColor)...
              +   A_loglike_BAYES_Beta(CChoice,CSafeOption,[x(:,1:2) x(:,9:10)],CNumDraws,CDraw,CDecisionColor)...
              +R_loglike2(RChoice,RProb,RSafeOption,x(:,1:2));
          


% 11 - MP in both tasks
ll{4} = @(x)    A_loglike_MPbeta_black(AChoice,AProb,ASafeOption,[x(:,1:2) x(:,5:8)],ANumDraws,ADraw,ADecisionColor)...
              +  A_loglike_MPbeta_black(CChoice,CProb,CSafeOption,[x(:,1:2) x(:,11:14)],CNumDraws,CDraw,CDecisionColor)...
              +R_loglike2(RChoice,RProb,RSafeOption,x(:,1:2));
     

theta0 = repmat([0.5 2 log(2) 0 0 0 0 0 log(2) 0 0 0 0 0],[N 1]);
 
    

%test the likelihood functions
tic
% Ambiguous task subjects
[ll{1}(theta0) ll{2}(theta0)]

toc


%% Some values that I keep referencing

numParams = size(theta0,2);
numModels = numel(ll);




%% set some simulation parameters

MH = 30; cMH = 0.1;
SimSize = 30000;
BurnIn =  10000;


ThetaM = theta0(1,:);


ThetaV = eye(numParams);

Theta = theta0;
MixProbs = ones(2)*0.5;
    % priors
    priorMix = [1 1 1 1]; % dirichlet prior, (ones==>multivariate uniform)
    Prior_MU = theta0(1,:);
    Prior_T  = 1;
    Prior_Sinv = eye(numParams);
    Prior_W = numParams+2;
 

 %% Arrays to dump stuff in to
 THETA    = zeros([N numParams SimSize]);
 THETA_MU = zeros([1 numParams SimSize]);
 THETA_V  = zeros([numParams numParams SimSize]);
 
 PRTYPE = zeros([N numModels SimSize]);
 ITYPE  = zeros([N 1 SimSize]);
 MIXPROBS = zeros([1 numModels SimSize]);
% URATE    = zeros([N 1 SimSize]);
 


%% Simulation

rand('state',42)
randn('state',4242)
mypool = parpool

tic

for ss = 1:SimSize
    % MH step for each type
    T = repmat(Theta,[1 1 numModels]);
    parfor tt = 1:numModels
        theta = T(:,:,tt);
        lt = ll{tt}(theta)+logmvnpdf(theta,ThetaM,ThetaV)';
    for mm = 1:MH
        prop = theta+cMH.*(chol(ThetaV)'*randn(numParams,N))';
       

    
        
        lp = ll{tt}(prop)+logmvnpdf(prop,ThetaM,ThetaV)';
        update = (lp-lt)>log(rand(N,1));
        theta(update==1,:) = prop(update==1,:);
        lt(update==1) = lp(update==1);
        % check for NaNs. Kill the program if it does
        if sum(isnan(lp))>0
            fprintf('ERROR: NaNs in the likelihoods. Banish them before believing the results\n')
            kill
        end
    
    end
    
        T(:,:,tt) = theta;
    end
    


   % evaluate types
    typePost = zeros(N,numModels);
    % log posterior probability of each type
    for tt = 1:numModels
        typePost(:,tt) = ll{tt}(T(:,:,tt))+log(MixProbs(tt));
    end
    % convert to actual probabilities
    typePost = typePost - repmat(max(typePost,[],2),[1 numModels]); % avoids exp(something large)
    typePost = exp(typePost)./repmat(sum(exp(typePost),2),[1 numModels]);
    % type categorical variable
    typeI = discretesampleJB(typePost);
    % update Theta to be equal to the value of the selected type;
    for tt = 1:4
        Theta(typeI==tt,:) = T(typeI==tt,:,tt);
    end
    % A NOTE ON DATA AUGMETATION HERE
    % Here we have a draw, Theta, from the posterior conditional (on everything
    % else). The MH process above makes sure that, say, for the Bayes Only type
    % (type 11), ALL parameters, not just the ones used in this model, are
    % drawn from the correct posterior conditional. Therefore, we don't have to
    % worry about further augmenting (in this case) the data with the
    % (irrelevant) MP parameters, because they are already drawn from the
    % correct distribution.

    % update mixing probabilities
    postMix = sum(repmat(typeI,[1 numModels])==repmat(1:numModels,[N 1]),1)+priorMix;
    mixProbGamma = random('Gamma',postMix,1);
    MixProbs = mixProbGamma./sum(mixProbGamma);

     % Hyperparametrers stage
         Post_MU = (Prior_T.*Prior_MU+N.*mean(Theta))./(Prior_T+N);
         Post_T = Prior_T+N;
         Post_W = N+Prior_W;
         S = zeros(numParams);
         for ii = 1:N
             S = S+(Theta(ii,:)-mean(Theta))'*(Theta(ii,:)-mean(Theta));
         end
         Post_S = Prior_Sinv + S + Prior_T.*N./(Prior_T+N).*(Prior_MU-mean(Theta))'*(Prior_MU-mean(Theta));
         % draw ThetaV ~ IW(Post_S,Post_W)
         ThetaV = eye(numParams)/wishrand(Post_W, Post_S);
         % Draw ThetaM ~ N(Post_MU, ThetaV/Post_W)
         ThetaM = Post_MU+(chol(ThetaV./Post_W)'*randn(size(ThetaM')))';%     ThetaM = mvnrnd(Post_MU,ThetaV./Post_W);

    % dump stuff into something to save

    THETA(:,:,ss)    = Theta;
    THETA_MU(:,:,ss) = ThetaM;
    THETA_V(:,:,ss)  = ThetaV;
    PRTYPE(:,:,ss)   = typePost;
    ITYPE(:,:,ss)    = typeI;
    MIXPROBS(:,:,ss) = MixProbs;
    
    disp([MixProbs])
    disp([int2str(ss) ' -- ' num2str(toc/ss*SimSize*(1-ss./SimSize)/3600) ' h to go'])
end
toc
%delete(mypool)
%% Store the simulation results (on a network drive accessable on JB's office computer)

fpath = ['C:\Users\jbland\Documents\LearningSims\' modelname]

save(fpath)






