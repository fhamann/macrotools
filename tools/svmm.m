function [Z,P]   = svmm(p,q,exprob,sigmaH,sigmaL,m,width,mH,mL)
% SVMM.M Stochastic volatility Markov Chain
%
% Usage:         [Z,P]   = svmm(p,q,exprob,sigmaH,sigmaL,m,width,mH,mL)
%
% Syntax: 
%   INPUTS
%      p       : probability of staying in low state
%      q       : probability of staying in high state
%      exprob  : exogenous probability matrix
%      sigmaH  : standard deviation of shocks in high state
%      sigmaL  : standard deviation of shocks in low state
%      m       : mean of the shocks
%      width   : "width" of discretized state space s   
%      mH      : mean of the shocks in high state
%      mL      : mean of the shocks in low state
%
%    OUTPUTS
%      P       : 4 by 4 transition probability for states 
%      Z       : 4 by 1 vector, nodes for Z
%


Prob = [p (1-p); (1-q) q];

P = kron(Prob,exprob);
N = length(P);


% sigmazH = sigmaH / sqrt(1-rho^2);
% sigmazL = sigmaL / sqrt(1-rho^2);

if nargin<7
width = 3;
end
if nargin<8
mH = m+width*sigmaH;
mL = m-width*sigmaL;
end


fil = sqrt((N/2-1))*sigmaL;
ZL  = linspace(-fil,fil,N/2)';
ZL  = ZL + mL;

fih = sqrt((N/2-1))*sigmaH;
ZH  = linspace(-fih,fih,N/2)';
ZH  = ZH + mH;

Z = [ZL ; ZH]'; 











