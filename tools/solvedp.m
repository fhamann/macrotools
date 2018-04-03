function [v,x,pstar] = solvedp(f,P,beta,alg,v)

% SOLVEDP  Solves stochastic dynamic program by state-space discretization
% Usage:
%         [v,x,pstar] = solvedp(f,P,beta,alg,v);
% or...
%         [v,x,pstar] = solvedp(f,P,beta);
% in the case for stochastic discount factor (when beta is n by m matrix)
%
% Syntax: Let m = # actions, n = # of all possible states, then
%   INPUTS
%      f       : n by m reward function
%      P       : mn by n transition probability matrix
%      alg     : algorithm used = 'policy' for policy function iteration
%                               = 'value'  for value fucntion iteration
%      v       : initial guess for value function (default is zero) 
%    OUTPUTS
%      v       : value function, v(s)
%      x       : optimal controls, x(s)
%      pstar   : optimal transition probability for states s, P(x(s))
%
% Adapted by Franz Hamann based on some old code by Fackler and Miranda.

% SET CONVERGENCE PARAMETER DEFAULTS
  maxit    = 5000;         % maximum number of iterations
  tol      = 10e-6;        % convergence tolerance
  prtiters = 50;           % print iterations every (prtiters>0) or not (0)

[n,m] = size(f);
if nargin<5 v=zeros(n,1);         % initial v [v=0]
else v=reshape(v,n,1);
end

% PERFORM POLICY OR VALUE FUNCTION ITERATIONS

if nargin<4
  if prtiters, 
      disp('Solving Bellman equation by value function iteration'); 
      %disp('Case stochastic discount factor');
  end
    for it=1:maxit
      vold = v;                             % store old value
      [v,x] = valmax(v,f,P,beta);           % update policy
      change = norm(v-vold);                % compute change
      if mod(it,prtiters)==0; fprintf ('%5i %10.1e\n',it,change); end
      if change<tol, break, end;            % convergence check
    end
    pstar = valpol(x,f,P,beta);
    if change>tol, warning('Failure to converge in solvedp'), end;
else
switch alg
  case 'policy'
    if prtiters, disp('Solving Bellman by policy function iteration'); end
    for it=1:maxit
       vold = v;                            % store old value
      [v,x] = valmax(v,f,P,beta);           % update policy  
      [pstar,fstar] = valpol(x,f,P,beta);   % induced P and f 
      v = (speye(n,n)-beta*pstar)\fstar;      % update value
      change = norm(v-vold);                % compute change
      if mod(it,prtiters)==0; fprintf ('%5i %10.1e\n',it,change); end
      if change<tol, break, end;            % convergence check
    end
    if change>tol, warning('Failure to converge in solvedp'), end;
  case 'value'
    if prtiters, disp('Solving Bellman by Value function iteration'); end
    for it=1:maxit
      vold = v;                             % store old value
      [v,x] = valmax(v,f,P,beta);           % update policy
      change = norm(v-vold);                % compute change
      if mod(it,prtiters)==0; fprintf ('%5i %10.1e\n',it,change); end
      if change<tol, break, end;            % convergence check
    end
    pstar = valpol(x,f,P,beta);
    if change>tol, warning('Failure to converge in solvedp'), end;
  otherwise
    error('algorithm must be policy or value')
  end
    
end
    
function [v,x] = valmax(v,f,P,beta)
[n,m]=size(f);
[v,x]=max(f+beta.*reshape(P*v,n,m),[],2);

function [pstar,fstar] = valpol(x,f,P,beta)
[n,m]=size(f); 
i=(1:n)';
fstar = f(n*(x-1)+i);
pstar = P(n*(x-1)+i,:);
