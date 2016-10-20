% RBC.M Solves standard real business cycle model with inelastic labor
%       supply. 
%
% Written by F. Hamann. Feel free to copy, change and distribute
 fprintf('\nReal business cycle model \n')

%% Parameters
 beta  =  0.98;                    % discount factor
 mu    =  2.0;                     % utility parameter
 alpha =  0.3;                     % production parameter
 delta =  0.1;                     % depreciation rate
 sigma =  0.35;                    % productivity shock volatility
 rho   =  0.5;                     % productivity shock persistence
  
%% Approximate z with nz discrete states AR(1)
 nz        =  9;                   % gridpoints in Z
 [z, prob] = rouwenhorst(nz,0,rho,sigma);
 z         = exp(z);
  
%% Construct state space
 nk       =  100;                  % gridpoints of K
 k        = linspace(0.1,20,nk);
 [Z,K,Kp] = gridmake(z',k',k');

 m     = nk;                       % number of actions
 n     = nk*nz;                    % number of states
 
%% Construct the return function, u
  C = Z.*K.^alpha + (1-delta)*K - Kp;
  u = C.^(1-mu)./(1-mu); u(C<=0) = NaN;
  u = reshape(u,n,m);
  
%% Construct the transition matrix
  P  = kron(speye(m,m),repmat(prob,m,1));
  
%% Solve the DP problem using policy iteration
  [v,x,pstar] = solvedp(u,P,beta,'policy'); clear P u C;
 
  d = ergdist(pstar);

  v = reshape(v,nz,nk);
  figure(1); surf(k,z,v);
  xlabel('k'); ylabel('z'); zlabel('v');
  title('Value function')

  d = reshape(d,nz,nk);
  figure(2); surf(k,z,d);
  xlabel('k'); ylabel('z'); zlabel('d');
  title('Stationary density')