function [ll] = A_loglike_BAYES_BetaRestricted(choices,safe,params,numdraws,draw,color)
%Returns the log-likelihood of Bayesian for compound
%choices. 
%Beta prior, truncated to [0.25,0.75]
[N,T] = size(choices);



draw(draw==-1)=0; % so Bayes' Rule works properly

gamma = params(:,1);
lambda = exp(params(:,2));

U33 = repmat(33.^(1-gamma)./(1-gamma),[1 T]);
U05 = repmat(5.^(1-gamma)./(1-gamma),[1 T]);
%Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));

% priors on the 3 urns'
N0 = repmat(exp(params(:,3)),[1 T]);
p0 = repmat(normcdf(params(:,4)),[1 T]);
% p0 is in terms of black, so if the decision is about white balls, adjust
% it. P0 below here refers to posterior probability of a win
p0(color==0) = 1-p0(color==0);


aBayes = draw+p0.*N0;
bBayes = numdraws-draw+(1-p0).*N0;

RejectPR = isnan(aBayes) | isnan(bBayes);

aBayes(RejectPR) = 1;
bBayes(RejectPR) = 1;

PrWin = (betainc(0.75,aBayes+1,bBayes)-betainc(0.25,aBayes+1,bBayes))...
           ./(betainc(0.75,aBayes,bBayes)-betainc(0.25,aBayes,bBayes))...
           .* (aBayes)./(aBayes+bBayes);

PrWin(isnan(PrWin)) = max(min(p0(isnan(PrWin)),0.75),0.25);


% EU of taking the risky option
EUrisky =sum(PrWin.*U33+(1-PrWin).*U05,3);
% U from taking the safe option
Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
% Difference in utility, Uchosen-Unotchosen
DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - EUrisky)./(U33-U05);

l = -log(1+exp(-DUchoices));

l(isnan(l)) = 0;
l(RejectPR) = -realmax;



ll = sum(l,2);


end