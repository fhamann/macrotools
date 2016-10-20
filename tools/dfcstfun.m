function [s_t,x_t] = dfcstfun(P,S,X,T,s,x0)

% DFCSTFUN Discrete forecast functions of a controled Markov Chain 
%
% Usage:
%          [s_t, x_t] = dfcstfun(P,S,X,x,T,s,x0)
% INPUTS:
%          P     : n by n transition probability matrix
%          S     : n by m matrix of n discrete states for m variables 
%          X     : n by k matrix of n discrete actions for k variables
%          x     : n by 1 policy function (optimal controls)
%          T     : number of simulated time periods
%          s     : 1 by m initial state *index* at t=0
%          x0    : 1 by k initial *values* of optimal controls at t=0 
%   where 
%          n     : total number of states in the Markov Chain
%          m     : total number of state variables
%          k     : total number of control variables
% OUTPUTS:
%          s_t   : T+1 by m path of states conditional on starting at s
%          x_t   : T+1 by k path of corresponding optimal controls 

 if sum(P')~=1; warning('Check accuracy of transition matrix'); end;

 [n,m] = size(S);   s_t = zeros(T,m);
 [~,k] = size(X);   x_t = zeros(T,k);
 
 f  = zeros(n,1); f(s)=1;        % initial degenerate probability at s 
 
 for t = 1:T;
  s_t(t,:) = f'*S;                % expected state at t
  x_t(t,:) = f'*X;                % expected control at t  
  f        = (f'*P)';             % update probabilities  
 end;
 
 x_t = [x0; x_t]; 
 s_t = [S(s,:); s_t];