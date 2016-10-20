function [P,f] = spmc(rho,rhoxy)

% SPMC.M  Simple persistence of 4-state 2-variable symmetric Markov Chain
% 
% Usage:
%                  [P,f] = spmc(rho,rhoxy)
%
%    INPUTS
%     rho      : autocorrelation of states (simple persistence)
%     rhoxy    : correlation coefficient between variables X and Y
%
%    OUTPUTS
%      P       : transition probability for states 
%      f       : ergodic distribution
%
% Explanation: there are two variables X and Y and one wants to create a
%              low-dimension symmetric Markov Chain: each variable has two 
%              states (high and low, for instance). In total, one would 
%              have four states {(x0,y0),(x0,y1),(x1,y0),(x1,y1)}. 
%              The correlation coefficient between X and Y is rhoxy.
%              In addition, the states of the process persist over time
%              obeying a "simple-persistence" rule as in Barton, David and
%              Fix (1962) "Persistence in a Chain of Multiple Events when
%              there is Simple Dependence", Biometrika Vol. 49 No. 3/4.

pi    = [(1+rhoxy)/4 0.5-(1+rhoxy)/4 0.5-(1+rhoxy)/4 (1+rhoxy)/4];
rho   = rho*ones(length(pi),1);
pi    = kron(1-rho,pi);
delta = diag(rho,0);
P     = pi+delta;
f     = markov(P);