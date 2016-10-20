function [spath,xpath,jpath] = imprespmarkov(P,S,x,T,s)

% IMPRESPMARKOV Impulse response of optimal controls of a Markov Chain 
%
% Usage:
%        [Spath,xpath,jpath] = imprespmarkov(P,S,x,T,s)
%
% INPUTS:
%          P     : n by n transition probability matrix
%          S     : n by m matrix of n discrete states for m variables  
%          x     : n by 1 policy function (optimal controls)
%          T     : number of simulated time periods
%          s     : initial state index
%
% where n is the total number of states in the Markov Chain
%       m is the total number of state variables
%
% OUTPUTS:
%          spath  : T by m path of states conditional on starting at s
%          xpath  : T by 1 path of corresponding optimal controls 
%          jpath  : T by 1 path of optimal state index

 if sum(P')~=1; 
     warning('Check transition probability matrix accuracy'); 
 end;

 [n,m] = size(S);
 spath = zeros(T,m);
 pi    = zeros(n,1); pi(s)=1;        % initial distribution
 
 for t = 1:T;
  spath(t,:) = pi'*S;                % expected state at t
  pi         = (pi'*P)';             % update the distribution pi  
 end;

 if nargout>1; xpath = x(getindex(spath,S)); end;  % optimal control at t
 if nargout>2; jpath = getindex(spath,S)   ; end;  % optimal state index
 