function [vs,xs,vc,xc,D] = solvedpos(fc,fs,prob,beta,phi,i0)

% SOLVEDP  Solves optimal stopping program by state-space discretization.
%
% Usage:
%
%   [vs,xs,vc,xc,D] = solvedpos(fc,fs,prob,beta,phi,i0)
%
% Syntax: Let m = # possible actions, n = # of possible states, then
%
%   INPUTS
%      fc     : n by m reward function if continues
%      fs     : n by m reward function if stops
%      prob   : transition probability matrix of exogenous states
%      beta   : discount factor
%      phi    : exogenous probability to re-entry from the outside
%      i0     : positions in X-grid for which debt is zero
%
%    OUTPUTS
%      vs      : value function if stops
%      xs      : policy function if stops
%      vc      : value function if continues
%      xc      : policy function if continues
%      D       : states in which vs>vc

% SET CONVERGENCE PARAMETER DEFAULTS

maxit    = 10000;            % maximum number of iterations
tol      = 10e-6;            % convergence tolerance, usually tol=sqrt(eps)
prtiters = 0;                % print iterations (1) or not (0)

[n,m] = size(fc);

vs = max(fs,[],2);
vc = max(fc,[],2);
v0 = getv0(i0,vc,n,m);
    
if prtiters, disp('\nSolving problem by value function iteration'); end
      
for it=1:maxit
    
    vsold=vs;    vcold=vc;
        
    Evc = kron(prob*(reshape(vc,m,n/m)'),ones(m,1));
    Evs = kron(prob*(reshape(vs,m,n/m)'),ones(m,1));
    Ev0 = kron(prob*(reshape(v0,m,n/m)'),ones(m,1));
      
    [vs,xs] = max(fs+phi*(beta.*Ev0)+(1-phi)*(beta.*Evs),[],2);
    [vc,xc] = max(fc+beta.*Evc,[],2);
      
    D = vs>vc  | isnan(vc)==1;   % Compute default decision
    vc(find(D))=vs(find(D));     % Max(vc,vs)
    
    v0  = getv0(i0,vc,n,m);     
    v = [vs vc];  vold = [vsold vcold]; if isnan(v), break, end;
      
    change = norm(v-vold);                
      
    if prtiters   
       fprintf ('%5i %10.1e\n',it,change)   
    end
    
    if change<tol, break, end;           
end

if change>tol,warning('No convergence'),else disp(' ... Solved!'); end;


function v0 = getv0(i0,vc,n,m)

v0 = reshape(vc,m,n/m);  
v0 = v0(i0,:);      
v0 = kron(v0,ones(m/size(v0,1),1));
v0 = reshape(v0,n,1); 