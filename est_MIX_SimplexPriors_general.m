%% simulates the posterior distribution assuming that:
% either MP of subjective bayes 
% parameters and type can be different between tasks
% Simplex priors used
% Parameters for simplex prior are the lower bounds of probability assigned
% to each urm
clear; clc;

load Exp1.mat

Experiment = 1;

doML = 0; % set ==1 to do individual ML estimates (slow and ultimately useless when there are a lot of parameters)

E1_DecisionColor_A= strcmp(E1_DecisionColor_A,'b');
E1_DecisionColor_R= strcmp(E1_DecisionColor_R,'b');
E1_DecisionColor_A2= strcmp(E1_DecisionColor_A2,'b');
E1_DecisionColor_C= strcmp(E1_DecisionColor_C,'b');


N = size(E1_Choice_C,1);

modelname = 'MixSimplexGeneral'

%% Set up the likelihood function


%parameter vector theta is:
%    parameter   task     type starting_value
% ------------------------------------------------
%01  gamma       all      all       0.5
%02  lambda      all      all       exp(2)
%03  p1          A        Bayes     0.5
%04  p2          A        Bayes     0.25
%05  alpha       A        MP        normcdf(1)
%06  pmin1       A        MP        normcdf(-2)   - bound urn Pr(1)
%07  pmin2       A        MP        normcdf(-2)   - bound urn Pr(2)
%08  pmin3       A        MP        normcdf(-2)   - bound urn Pr(3)
%09  p1          C        Bayes     0.5
%10  p2          C        Bayes     0.25
%11  alpha       C        MP        0.7
%12  pmin1       C        MP        normcdf(-2)   - bound urn Pr(1)
%13  pmin2       C        MP        normcdf(-2)   - bound urn Pr(2)
%14  pmin3       C        MP        normcdf(-2)   - bound urn Pr(3)


theta0 = repmat([0.5 2 0 0 1 -2 -2 -2 0 0 1 -2 -2 -2],[N 1]);

l= A_loglike_BAYES_Simplex(E1_Choice_C,E1_SafeOption_C,theta0(:,1:4),E1_NumDraws_C,E1_Draw_C,E1_DecisionColor_C);
l= A_loglike_MP_Simplex3([E1_Choice_A E1_Choice_A2],[E1_SafeOption_A E1_SafeOption_A2],[theta0(:,1:2) theta0(:,5:8)],[E1_NumDraws_A E1_NumDraws_A2],[E1_Draw_A E1_Draw_A2],[E1_DecisionColor_A E1_DecisionColor_A2]);

% 00 - Bayes in both tasks
ll{1} = @(x)    A_loglike_BAYES_Simplex([E1_Choice_A E1_Choice_A2],[E1_SafeOption_A E1_SafeOption_A2],x(:,1:4),[E1_NumDraws_A E1_NumDraws_A2],[E1_Draw_A E1_Draw_A2],[E1_DecisionColor_A E1_DecisionColor_A2])...
              +A_loglike_BAYES_Simplex(E1_Choice_C,E1_SafeOption_C,[x(:,1:2) x(:,9:10)],E1_NumDraws_C,E1_Draw_C,E1_DecisionColor_C)...
              +R_loglike2(       E1_Choice_R,E1_Prob_R,E1_SafeOption_R,x(:,1:2));
% 01 - A-Bayes & C-MP
ll{2} = @(x)    A_loglike_BAYES_Simplex([E1_Choice_A E1_Choice_A2],[E1_SafeOption_A E1_SafeOption_A2],x(:,1:4),[E1_NumDraws_A E1_NumDraws_A2],[E1_Draw_A E1_Draw_A2],[E1_DecisionColor_A E1_DecisionColor_A2])...
              +A_loglike_MP_Simplex3(E1_Choice_C,E1_SafeOption_C,[x(:,1:2) x(:,11:14)],E1_NumDraws_C,E1_Draw_C,E1_DecisionColor_C)...
              +R_loglike2(       E1_Choice_R,E1_Prob_R,E1_SafeOption_R,x(:,1:2));

% 10 - A-MP & C-Bayes
ll{3} = @(x)    A_loglike_MP_Simplex3([E1_Choice_A E1_Choice_A2],[E1_SafeOption_A E1_SafeOption_A2],[x(:,1:2) x(:,5:8)],[E1_NumDraws_A E1_NumDraws_A2],[E1_Draw_A E1_Draw_A2],[E1_DecisionColor_A E1_DecisionColor_A2])...
              +   A_loglike_BAYES_Simplex(E1_Choice_C,E1_SafeOption_C,[x(:,1:2) x(:,7:8)],E1_NumDraws_C,E1_Draw_C,E1_DecisionColor_C)...
              +R_loglike2(       E1_Choice_R,E1_Prob_R,E1_SafeOption_R,x(:,1:2));
          


% 11 - MP in both tasks
ll{4} = @(x)    A_loglike_MP_Simplex3([E1_Choice_A E1_Choice_A2],[E1_SafeOption_A E1_SafeOption_A2],[x(:,1:2) x(:,5:8)],[E1_NumDraws_A E1_NumDraws_A2],[E1_Draw_A E1_Draw_A2],[E1_DecisionColor_A E1_DecisionColor_A2])...
              +   A_loglike_MP_Simplex3(E1_Choice_C,E1_SafeOption_C,[x(:,1:2) x(:,11:14)],E1_NumDraws_C,E1_Draw_C,E1_DecisionColor_C)...
              +R_loglike2(       E1_Choice_R,E1_Prob_R,E1_SafeOption_R,x(:,1:2));

%test the likelihood functions
tic
for cc = 1:4;
    ll{cc}(theta0)
end
toc



%% A big note for elegant coding:
% Everything below here is the same for the mixture models. It should be in
% the same file, and it will be on the next iteration. Future James: you
% are on notice!

%% Some values that I keep referencing

numParams = size(theta0,2);
numModels = numel(ll);




%% set some simulation parameters

MH = 30; cMH = 0.1;
SimSize = 30000;
BurnIn =  10000;


% starting values
ThetaM = theta0(1,:);                       % mean parameter vector
ThetaV = eye(numParams);                    % variance of parameter vector
MixProbs = ones(1,numModels)./numModels;    % mixing probabilities

% priors
priorMix = ones(1,numModels); % dirichlet prior, (ones==>multivariate uniform)
 Prior_MU = theta0(1,:);
 Prior_T  = 1;
 Prior_Sinv = eye(numParams);
 Prior_W = numParams+2;
 

%% starting values

Theta = theta0;

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
%mypool = parpool

tic

for ss = 1:SimSize
    % MH step for each type
    T = repmat(Theta,[1 1 numModels]);
    for mm = 1:MH
        prop = T + repmat(cMH.*(chol(ThetaV)'*randn(size(Theta')))',[1 1 numModels]);
        %lt = zeros(N,numModels);
        %lp = zeros(N,numModels);
        update = zeros(N,1,numModels);
    for tt = 1:numModels; % CAN I PARALELLIZE THIS STEP? it is the slowest!!!
        lt = ll{tt}(T(:,:,tt))+logmvnpdf(T(:,:,tt),ThetaM,ThetaV)';
        lp = ll{tt}(prop(:,:,tt))+logmvnpdf(prop(:,:,tt),ThetaM,ThetaV)';
        
        %end
        %update = (lp-lt)>log(rand(N,numModels));
        %for tt = 1:numModels
        %    T(update(:,tt),:,tt)= prop(update(:,tt),:,tt);
        %end
        update(:,1,tt) = (lp-lt)>log(rand(N,1));
        T(update(:,1,tt)==1,:,tt) = prop(update(:,1,tt)==1,:,tt);
        
        % check for NaNs. Kill the program if it does
        if sum(isnan(lp))>0
            fprintf('ERROR: NaNs in the likelihoods. Banish them before believing the results\n')
            kill
        end
    
    end
    
        %T(repmat(update==1,[1 numParams 1])) = prop(repmat(update==1,[1 numParams 1]));
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
    %  plot(permute(MixProbs,[3 2 1]))
end
toc
%delete(mypool)
%% Store the simulation results (on a network drive accessable on JB's office computer)

fpath = ['C:\Users\jbland\Documents\LearningSims\' modelname]

save(fpath)






