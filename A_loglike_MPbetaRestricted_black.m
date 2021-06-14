function [ll] = A_loglike_MPbetaRestricted_black(choices,prob,safe,params,numdraws,draw,color)
%Returns the log-likelihood of the multiple priors model for ambiguous
%choices. Set of Beta priors
% framed in terms of black draws


% later, when assessing the worst-case scenario, for white decisions I need
% to find the maximum Pr(black|s,mu), not the minimum.

% Restriction: since in this model the subject knows that p\in[0.25,0.75]
% they have priors of the form p ~ truncatedBeta(p0N0,(1-p0)N0,[0.25,0.75])

gs = 21; % size of gridsearch between p0L and p0R

draw(draw==-1)=0; % so Bayes' Rule works properly

% A note on probabilities: I (JB) define p0 in terms of Pr(black).
% Therefore I need to change the variable "draw" to 1-draw when the
% decision color is white:

draw(color==0) = numdraws(color==0)-draw(color==0);

black = draw; black(color==0)=numdraws(color==0)-draw(color==0);
black(isnan(black))=0; numdraws(isnan(numdraws))=0;
[N,T] = size(choices);


gamma = params(:,1);
lambda = exp(params(:,2));
a = repmat(normcdf(params(:,3)),[1 size(prob,2)]); % alpha parameter
N0 = exp(params(:,4));
% For experiment 2, for some estimations we impose the restriction that p0R=1-p0L
% We let the function know this by putting NaNs in the 6th column of
% parameters.
if sum(isnan(params(:,6)))>0
    p0L = 0.5.*normcdf(params(:,5));
    p0R = 1-p0L;
else
    pmid = normcdf(params(:,5));
    p0L = pmid-pmid.*normcdf(params(:,6));
    p0R = pmid + (1-pmid).*normcdf(params(:,6));
end



U33 = repmat(33.^(1-gamma)./(1-gamma),[1 T]);
U05 = repmat(5.^(1-gamma)./(1-gamma),[1 T]);
Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));


p0grid = repmat(permute(linspace(0,1,gs),[1 3 2]),[N T 1]);
p0grid = repmat(p0L,[1 T gs]) + repmat(p0R-p0L,[1 T gs]).*p0grid;

BinCoef = repmat(factorial(numdraws)./(factorial(black).*factorial(numdraws-black)),[1 1 gs]); % Binomial coefficient


aX = repmat(black,[1 1 gs]) +repmat(N0,[1 T gs]).*p0grid;
bX = repmat(numdraws-black,[1 1 gs]) + repmat(N0,[1 T gs]).*(1-p0grid);

% When N0 is large enough for dirac priors
InfBeta = isinf(aX) | isinf(bX);
aX(InfBeta) = 1;
bX(InfBeta) = 1;


% restriction goes in here - truncated beta distribution
priorLike =   BinCoef ...
    .* beta(aX,bX)...
    .*(betainc(0.75,aX,bX)-betainc(0.25,aX,bX));
priorLikeSub = BinCoef.*p0grid.^repmat(black,[1 1 gs]).*(1-p0grid).^repmat(numdraws-black,[1 1 gs]);
priorLike(InfBeta) = priorLikeSub(InfBeta);
 
maxL = max(priorLike,[],3);

dropPrior = priorLike <= repmat((a.*maxL),[1 1 gs]); % priors that are not discarded



maxp0kept = max(p0grid - dropPrior.*realmax,[],3);   
minp0kept = min(p0grid + dropPrior.*realmax,[],3);
% here we run into a problem if N0 is large: we drop all of the priors. In
% this case the subject is effectively using dirac priors, and for
% alpha\neq0 will drop all prior except those with p0 = black/numdraws.
% Therefore the only remaining prior here is either black/numdraws, p0L, or
% p0R.
pNlarge = max(min(black./numdraws,repmat(p0R,[1 T])),repmat(p0L,[1 T]));

    maxp0kept = maxp0kept.*(numdraws>0)+max(p0grid,[],3).*(1-(numdraws>0));
    minp0kept = minp0kept.*(numdraws>0)+min(p0grid,[],3).*(1-(numdraws>0));


maxp0kept(maxp0kept<0) = pNlarge(maxp0kept<0);
minp0kept(minp0kept>1) = pNlarge(minp0kept>1);
    
%Ep33black = (black+minp0kept.*repmat(N0,[1 T]))./(numdraws+repmat(N0,[1 T]));


aBlack = black+minp0kept.*repmat(N0,[1 T]);
bBlack = numdraws-black+(1-minp0kept).*repmat(N0,[1 T]);

RejectPR = isnan(aBlack) | isnan(bBlack);

aBlack(RejectPR) = 1;
bBlack(RejectPR) = 1;

Ep33black = (betainc(0.75,aBlack+1,bBlack)-betainc(0.25,aBlack+1,bBlack))...
           ./(betainc(0.75,aBlack,bBlack)-betainc(0.25,aBlack,bBlack))...
           .* (aBlack)./(aBlack+bBlack);
Ep33black(isnan(Ep33black)) = max(min(minp0kept(isnan(Ep33black)),0.75),0.25);
       
aWhite =  numdraws-black+(1-minp0kept).*repmat(N0,[1 T]);
bWhite = black+minp0kept.*repmat(N0,[1 T]);     

Ep33white = (betainc(0.75,aWhite+1,bWhite)-betainc(0.25,aWhite+1,bWhite))...
           ./(betainc(0.75,aWhite,bWhite)-betainc(0.25,aWhite,bWhite))...
           .* (aWhite)./(aWhite+bWhite);
Ep33white(isnan(Ep33white)) = max(min(1-minp0kept(isnan(Ep33white)),0.75),0.25);

PrWin = color.*Ep33black + (1-color).*Ep33white;


EUrisky =PrWin.*U33+(1-PrWin).*U05;
Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));

DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - EUrisky)./(U33-U05);

l = -log(1+exp(-DUchoices));

l(isnan(l)) = 0;
l(RejectPR) = -realmax;

ll = sum(l,2);


% % some tests to comment out when running the simulation
% disp(mean(dropPrior(:)))
% 
% figure;
% hold all
%     plot((DUchoices(draw>0)),PrWin(draw>0),'xk')
%     
% hold off


end

