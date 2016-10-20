function d = ergdist(P,d)

% ERGDIST.M Ergodic distribution of a n-state Markov Chain
%
% Usage:            d = ergdist(P,d)
%
% Input:
%          P  : n by n transition probability matrix
%          d  : n by 1 initial guess of invariant distribution
%
% where n is the total number of states in the Markov Chain.
%
% Ouput:
%          d  : n by 1 ergodic distribution
%
% Note: Does not work for Markov Chains with transient states.
%       Follows Ljundqvist and Sargent (2000). 

% if sum(P')~=1; warning('Check accuracy of transition matrix'); end;

n = length(P);

if nargin<2
  d = (1/n)*ones(n,1);  % Start with uniform dist., if no seed provided.
end

maxit = 500;
tol   = sqrt(eps);

for it = 1:maxit
  dnext  = (d'*P)';
  change = max(abs(dnext-d));
  %fprintf ('%5i %10.1e\n',it,change) 
  if change<tol; break; end;
  d  = dnext;
end  
  
