% CHRIST_FISHER Solves Real Business Cycle Model (inelastic labor supply)
%               with Reversible and Irreversible Investment
%
% See Christiano and Fisher (1994) "Algorithms for Solving Dynamic Models
% With Occasionally Binding Constraints" FRBM, Staff Report 171.
%
% Written by F. Hamann. Feel free to copy, change and distribute.

% Model parameters
  irrev =  1;                                     % irreversible investment: 1 (yes), 0 (no)
  beta  =  1.03^(-0.25);                          % discount factor
  mu    =  1.0;                                   % utility parameter
  alpha =  0.3;                                   % production parameter
  delta =  0.02;                                  % depreciation rate
  sigma =  0.22;                                  % production shock volatility
  rho   =  0;                                     % production shock persistence
  n1 =  1000;                                     % number of endogenous states
  n2 =  2;                                        % number of exogenous states
  m  =  n1;                                       % number of actions
  n  = n1*n2;                                     % total number of states

% Approximate productivity shocks with n2 states Markov chain
[a, prob] = rouwenhorst(n2,0,rho,sigma)
 a = exp(a);

% Construct state space
  kmin  = 16.9;
  kmax  = 55.1;
  k = linspace(kmin,kmax,n1);
  S = gridmake(a',k');
  [A,K,Knext] = gridmake(a',k',k');
  % Here is where the constraint is imposed
  if irrev == 1; disp('Model with Irreversible Investment')
    i = find(Knext < (1-delta)*K);
    Knext(i) = (1-delta)*K(i);
  else disp('Model with Reversible Investment')
  end

% Construct the reward function, u
  c = A.*K.^alpha + (1-delta)*K - Knext;
  if mu ==1;
    u = log(c);
  else
    u = c.^(1-mu)./(1-mu);
  end;

  j = find(c<0); % Consumption cannot be negative
  u(j) = NaN;
  u = reshape(u,n,m);
  
clear c i j A K Knext; 

% Construct the transition matrix
P = kron(speye(n1,n1),repmat(prob,n1,1));

% Solve infinite-horizon model via 'policy' or 'value' function iteration
  [v,x,Popt] = solvedp(u,P,beta,'policy');

clear u P; 

% Computes ergodic distribution
   piss = ergdist(Popt);

% Graphics
% First: reshaping for plots
  kopt    = k(x);
  vr      = reshape(v,n2,n1)';
  koptr   = reshape(kopt,n2,n1)';
  invest  = kopt'-(1-delta)*S(:,2);
  investr = reshape(invest,n2,n1)';
% Then the plots:
  figure(1); plot(k,vr); title('Value Function')
  xlabel('k'); ylabel('v');

  figure
  figure(2); plot(S(:,2),piss); title('Ergodic Distribution of K')
  xlabel('k'); ylabel('probability')

  figure
  figure(3); plot(k,investr); title('Investment Policy')
  xlabel('k'); ylabel('investment')
