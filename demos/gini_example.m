% Gini example for income distribution

addpath ~/Dropbox/matlab/inequality_package

nz     = 9;
zmean  = 1;
rhoz   = 0.3;                             
sigmaz = 0.75;                           

% Create exogenous transition matrix
[z, P] = rouwenhorst(nz,zmean,rhoz,sigmaz)
pzs    = markov(P);
y      = exp(z')

p = cumsum(pzs)

g = ginicoeff(p, y)
h = lorenzcurve(p,y);