function [ll] = R_loglike2(choices,prob,safe,params)
% inputs:
%   choices: subjects' binary choices (1=safe)
%   prob:    probability of $33 outcome
%   safe:    payoff of safe choice



gamma = params(:,1);
lambda = exp(params(:,2));

[N,T] = size(choices);


% 2016-10-03 something is going wrong here. Check later
Usafe = safe.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
U33   = 33.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
U05   = 5.^(1-repmat(gamma,[1 T]))./(1-repmat(gamma,[1 T]));
DUchoices = (1-2.*choices).*repmat(lambda,[1 T]).*(Usafe - prob.*U33-(1-prob).*U05)./(U33-U05);

l = -log(1+exp(-DUchoices));

l(isnan(choices)) = 0;

ll = sum(l,2);






end

