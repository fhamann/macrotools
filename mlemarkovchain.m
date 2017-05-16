function P = mlemarkovchain(X)

% MLEMARKOVCHAIN Maximum Likelihood Estimation for Markov Chains
%
% Usage:     P = mlemarkovchain(X)
%
% Inputs:
%        X: T by 1 vector of realizations of m-state MC, X={1,2,...,m}
%
% Outputs:
%       P: m by m MLE of MC computed as P(i,j) = N(i,j)/sum_j(N(i,j))
%          where N(i,j) is the number of times i is followed by j in X
%          and sum_j denotes the sum over j.


m = max(X);
n = length(X)-1;
P = zeros(m,m);

for t = 1:n
  P(X(t),X(t+1)) = P(X(t),X(t+1)) + 1;
end

for i = 1:m
  P(i,:) = P(i,:)/sum(P(i,:));
end