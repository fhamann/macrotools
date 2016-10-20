% Gini example for income distribution

addpath /Users/franz/Dropbox/matlab/inequality_package

nz = 9;
zmean  = 1;
rhoz   = 0.3;                             
sigmaz = 0.75;                           

% Create exogenous transition matrix
[z, prob] = rouwenhorst(nz,zmean,rhoz,sigmaz)
pzs       = markov(prob);
y         = exp(z');

p = cumsum(pzs)

g = ginicoeff(p, y)
h = lorenzcurve(p,y);