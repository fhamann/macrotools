function P = otpm(x,prob,ns1,ns2,nz)

% OTPM Optimal Transition Probability Matrix
%
% Usage:
%           P = otpm(x,prob,nb,ns,nz)
%
%   INPUTS
%       x      : optimal policy function (nb*ns*nz by 1  vector)
%       prob   : transition probability matrix for exogenous process
%       nb     : # of grid pts in the first endogenous state variable 
%       ns     : # of grid pts in the second endogenous state variable 
%       nz     : # of grid pts in the space of exogenous state variables 
%
%   OUTPUT
%       P      : Optimal transtion matrix P(s,z)


n   = length(x);
PP  = [prob;zeros((ns1*ns2-1)*nz,nz)];
PP2 = reshape(PP,nz,n);
PP3 = kron(PP2,ones(ns1*ns2,1));

P = zeros(n);

for i=1:n
    P(i,x(i):n) = PP3(i,1:n-x(i)+1);
end