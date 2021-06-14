function [ll] = A_loglike_MPbeta_black(choices,prob,safe,params,numdraws,draw,color)
%Returns the log-likelihood of the multiple priors model for ambiguous
%choices. Set of Beta priors
% framed in terms of black draws


% later, when assessing the worst-case scenario, for white decisions I need
% to find the maximum Pr(black|s,mu), not the minimum.


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


priorLike =   BinCoef.*beta(repmat(black,[1 1 gs]) +repmat(N0,[1 T gs]).*p0grid, repmat(numdraws-black,[1 1 gs]) + repmat(N0,[1 T gs]).*(1-p0grid)) ...
                     ./ beta(                     repmat(N0,[1 T gs]).*p0grid,                                   repmat(N0,[1 T gs]).*(1-p0grid));


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
    
Ep33black = (black+minp0kept.*repmat(N0,[1 T]))./(numdraws+repmat(N0,[1 T]));
Ep33white = (numdraws-black+(1-maxp0kept).*repmat(N0,[1 T]))./(numdraws+repmat(N0,[1 T]));
PrWin = color.*Ep33black + (1-color).*Ep33white;


EUrisky =PrWin.*U33+(1-PrWin).*U05;
Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));

DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - EUrisky)./(U33-U05);

l = -log(1+exp(-DUchoices));

l(isnan(l) | isnan(choices)) = 0;


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

