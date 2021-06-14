function [ll] = A_loglike_MP_Simplex3(choices,safe,params,numdraws,draw,color)
%Returns the log-likelihood of Bayesian for compound
%choices. 
%Beta prior
[N,T] = size(choices);



draw(draw==-1)=0; % so Bayes' Rule works properly


black = draw.*color + (numdraws-draw).*(1-color);

gamma = params(:,1);
lambda = exp(params(:,2));
alpha = normcdf(params(:,3));
pmin2  = normcdf(params(:,4));
pmin13 = (1-pmin2).*normcdf(params(:,5));
pmin1  = pmin13.*normcdf(params(:,6));
pmin3  = pmin13.*(1-normcdf(params(:,6)));

% Likelihood of priors at corners of the simplex
% 3rd dimension: corner of simplex with most probability.
% i.e. 2nd slice corresponds to (pmin,1-2pmin,pmin)
MinProbs = zeros(N,T,3);% repmat(pmin,[1 T 3]);
MinProbs(:,:,1) = repmat(pmin1,[1 T]);
MinProbs(:,:,2) = repmat(pmin2,[1 T]);
MinProbs(:,:,3) = repmat(pmin3,[1 T]);
CornerProbs3 = zeros(N,T,3,3);
LCorner = zeros(N,T,3);
for kk = 1:3
    CornerProbs = MinProbs;
    CornerProbs(:,:,kk) = 1-sum(MinProbs(:,:,(1:3)~=kk),3);
    CornerProbs3(:,:,:,kk) = CornerProbs;
    LCorner(:,:,kk) = sum(CornerProbs.*(repmat(permute(1:3,[1 3 2])./4,[N T 1])).^repmat(black,[1 1 3]).*(1-repmat(permute(1:3,[1 3 2])./4,[N T 1])).^repmat(numdraws-black,[1 1 3]),3);
end
% L* and the corner prior that maximizes it
[Lmax,maxK] = max(LCorner,[],3);
% most likely prior
maxP = zeros(N,T,3);
for cc = 1:3
    maxP = maxP+repmat((maxK==cc),[1 1 3]).*CornerProbs3(:,:,:,cc);
end

% identify corners of remaining set of priors
% That is, for each corner, since the likelihood is a linear function of
% the prior, the corner solutions of mu to L(mu) = alpha L* can be solved as:
% K(Lmax-LCorner) = (alpha*Lmax-Lcorner)
% i.e. alpha*Lmax = LCorner + K(Lmax-LCorner)
% Then: PCornerNew = PCorner + K (PMax-PCorner) = K*PMax+(1-K)*PCorner
K = (repmat(alpha,[1 T 3]).*repmat(Lmax,[1 1 3])-LCorner)./(repmat(Lmax,[1 1 3])-LCorner);
% If we get a K inside the unit interval, then we update this corner of the
% set of priors
updateK = K>0 & K<1;
updateK(K==-inf)=0;


updateK = repmat(permute(updateK,[1 2 4 3]),[1 1 3 1]);
PCornerInt = zeros(size(CornerProbs3));
for cc = 1:3
    PCornerInt(:,:,:,cc) = repmat(K(:,:,cc),[1 1 3]).*maxP + repmat(1-K(:,:,cc),[1 1 3]).*CornerProbs3(:,:,:,cc);
end

CornerProbs3(updateK) = PCornerInt(updateK);

% at this point, CornerProbs3 contains the priors at the corners of the set
% of priors. Since the problem is *very linear*, we only need to evaluate
% these priors to determine the subject's objective function.

 
 BinCoef = factorial(numdraws)./factorial(numdraws-black)./factorial(black);
 
 LikeUrn =  repmat(BinCoef,[1 1 3]).*repmat(permute(1:3,[1 3 2])./4,[N T 1]).^repmat(black,[1 1 3]).*(1- repmat(permute(1:3,[1 3 2])./4,[N T 1])).^repmat(numdraws-black,[1 1 3]);
 
 
 
 PostPrBlack = zeros(N,T,3);
 for kk = 1:3
     PostPrBlack(:,:,kk) = sum(repmat(permute(1:3,[1 3 2])./4,[N T 1]).*(CornerProbs3(:,:,:,kk).*LikeUrn)./repmat(sum(CornerProbs3(:,:,:,kk).*LikeUrn,3),[1 1 3]),3);
 end
 
 PostPrWin = PostPrBlack;
 PostPrWin(repmat(color,[1 1 3])==0) = 1-PostPrBlack(repmat(color,[1 1 3])==0);
 
 WorstPr = min(PostPrWin,[],3);
 
 U33 = repmat(33.^(1-gamma)./(1-gamma),[1 T]);
 U05 = repmat(5.^(1-gamma)./(1-gamma),[1 T]);
 EUrisky =sum(WorstPr.*U33+(1-WorstPr).*U05,3);
 Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
 DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - EUrisky)./(U33-U05);
 l = -log(1+exp(-DUchoices));
 ll = sum(l,2);
 
 
 
 
 


end