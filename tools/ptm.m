function pi = ptm(rho,rhoen);

% Creates 2 point Markov Chain shocks and PTM
% Usage:
%          pi = ptm(rho,rhoen)
%
% rho   : autocorrelation of shocks
% rhoen : correlation between shocks

pi    = [(1+rhoen)/4 0.5-(1+rhoen)/4 0.5-(1+rhoen)/4 (1+rhoen)/4];
rho   = rho*ones(4,1);
pi    = kron(1-rho,pi);
delta = diag(rho,0);
pi    = pi+delta;