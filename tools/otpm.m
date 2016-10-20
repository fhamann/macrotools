function P = otpm(x,prob,nb,ns,nz)

% OTPM Optimal Transition Probability Matrix
%
% Usage:
%           P = otpm(x,prob,nb,ns,nz)
%
%   INPUTS
%       x      : optimal policy function (column vector)
%       prob   : transition probability matrix for exogenous process
%       nb     : Number of grid points in the first endogenous variable 
%       ns     : Number of grid points in the first endogenous variable 
%       nz     : Number of grid points in the space of exogenous variables 
%
%   OUTPUT
%       P      : Optimal transtion matrix


n   = length(x);
PP  = [prob;zeros((nb*ns-1)*nz,nz)];
PP2 = reshape(PP,nz,n);
PP3 = kron(PP2,ones(nb*ns,1));

P = zeros(n);

for i=1:n
    P(i,x(i):n) = PP3(i,1:n-x(i)+1);
end