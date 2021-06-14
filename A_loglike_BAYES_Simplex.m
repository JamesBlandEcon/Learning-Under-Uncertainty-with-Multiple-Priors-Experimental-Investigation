function [ll] = A_loglike_BAYES_Simplex(choices,safe,params,numdraws,draw,color)
%Returns the log-likelihood of Bayesian for compound
%choices. 
%Simplex prior
[N,T] = size(choices);



draw(draw==-1)=0; % so Bayes' Rule works properly


black = draw.*color + (numdraws-draw).*(1-color);

gamma = params(:,1);
lambda = exp(params(:,2));

% prior: p(:,x) = probability assigned to there being x black balls in the urn
p = zeros(N,3);
p(:,2)     = normcdf(params(:,4));
p(:,1)     = (1-p(:,2)).*normcdf(params(:,3));
p(:,3)     = 1-p(:,1)-p(:,2);

% posterior distribution of urns
q = zeros(N,T,3);
for pp = 1:3
    q(:,:,pp) = repmat(pp./4,[N T]).^black.*repmat(1-pp./4,[N T]).^(numdraws-black).*repmat(p(:,pp),[1 T]);
end
q = q./repmat(sum(q,3),[1 1 3]);

% posterior expectectation of Pr(Black)
EpBlack = sum(q.*repmat(permute((1:3)./4,[3 1 2]),[N T]),3);

% Probability of winning the risky bet
PrWin = EpBlack.*color + (1-EpBlack).*(1-color);

U33 = repmat(33.^(1-gamma)./(1-gamma),[1 T]);
U05 = repmat(5.^(1-gamma)./(1-gamma),[1 T]);
EUrisky =sum(PrWin.*U33+(1-PrWin).*U05,3);
Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - EUrisky)./(U33-U05);
l = -log(1+exp(-DUchoices));
ll = sum(l,2);


% U33 = repmat(33.^(1-gamma)./(1-gamma),[1 T]);
% U05 = repmat(5.^(1-gamma)./(1-gamma),[1 T]);
% Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
% 
% % priors on the 3 urns
% 
% 
% 
% 
% p0(color==0) = 1-p0(color==0);
% 
% PrWin = (draw+p0.*N0)./(numdraws+N0); %Posterior expectation of winning probability
% 
% % EU of taking the risky option
% EUrisky =sum(PrWin.*U33+(1-PrWin).*U05,3);
% % U from taking the safe option
% 
% % Difference in utility, Uchosen-Unotchosen
% DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - EUrisky)./(U33-U05);
% 
% l = -log(1+exp(-DUchoices));
% 
% 
% 
% ll = sum(l,2);


end