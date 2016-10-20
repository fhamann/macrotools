function Ps  = betabinomial(s,n)

% BETABINOMIAL  Beta-binomial model for a 2-state Markov chain 
%
% Usage:      Ps  = betabinomial(s,n)
%
% Input:
%
%       s  : T by 1 column vector of positions of the states
%       n  : 2 by 2 matrix of initial counters 
%
% Output:
%
%       Ps : sequence t=1,...,T of transition probability matrices as a
%            function of the vector s conditional on past history up to t


T  = length(s);
ns = length(n);

Ps = zeros(ns,ns,T);

for t=2:T   
    for i=1:2
       for j=1:2
           if (s(t) == j) && (s(t-1)==i)
               n(i,j,t) = n(i,j,t-1)+1;
           else
               n(i,j,t) = n(i,j,t-1);
           end
       end
    end
end

nhh = squeeze(n(1,1,:));
nhl = squeeze(n(1,2,:)); 
nlh = squeeze(n(2,1,:));
nll = squeeze(n(2,2,:));

E_phh = nhh./(nhh+nhl);
E_phl = nhl./(nhh+nhl);
E_pll = nll./(nll+nlh);
E_plh = nlh./(nll+nlh);

for t=1:T
  Ps(:,:,t)=[E_pll(t) 1-E_pll(t) ; 1-E_phh(t) E_phh(t)];
end