function [Z,P]   = svmc(p,q,pL,pH,e,m,a)

% SVMC.M Stochastic volatility 4-states Markov Chain
%
% Usage:            P   = svmc(p,q,pL,pH)
%
% Syntax: 
%   INPUTS
%      p       : probability of staying in low output
%      q       : probability of staying in high output
%      pL      : probability of staying in low volatility regime
%      pH      : probability of staying in high volatility regime 
%      e       : scalar, size of epsilons
%      m       : scalar, value of mean for random variable
%      a       : scalar, value of excess volatility, a>1
%
%    OUTPUTS
%      P       : 4 by 4 transition probability for states 
%      Z       : 4 by 1 vector, nodes for Z
%
% Explanation: there are two regimes high volatility and low volatility.
%              One wants to create a low-dimension Markov Chain for output
%              in the two regimes: in each regime, output can be high or
%              low. Every period, if output is low, next period will be low
%              low with probability q. On the contrary, if it is high the
%              probability of staying high next period is p. 
%              Output dispersion in the high volatility regime varies
%              between {yH,YH}, while output dispersion in the low
%              volatility regime varies between {yL,YL}. In total, one ends 
%              up with four states {yH < yL < YL < YH}. The probability
%              that the system stays in the low volatility regime is pL,
%              while the probability of staying in the high volatility 
%              regime is pH. This code simulates this type of environment.
%              We label this as "stochastic volatility". 

 P  = [   pH*p       (1-pH)*p     (1-pH)*(1-p)   pH*(1-p)   ;
        (1-pL)*p       pL*p        pL*(1-p)    (1-pL)*(1-p) ;
       (1-pL)*(1-q)   pL*(1-q)       pL*q       (1-pL)*q    ;
         pH*(1-q)   (1-pH)*(1-q)   (1-pH)*q       pH*q      ];

 N  = length(P);
 fi = sqrt(N-1)*e;
 Z  = linspace(-fi,fi,N)';
 Z  = Z + m;
 
 if nargin>6; Z(1) = -a*fi; Z(end) = a*fi; end