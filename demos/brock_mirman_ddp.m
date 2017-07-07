% BROCK_MIRMAN_DDP Brock and Mirman (1972) optimal stochastic growth model
%                  by discrete dynamic programming
%
% Written by F. Hamann. Feel free to copy, change and distribute
 fprintf('\nBrock and Mirman model \n')

%% Parameters
 gamma	= 1;        % risk aversion 
 beta	= 0.98;     % discount factor
 alpha  = 0.3;      % capital share

%% Markov chain for z
 ne = 8;
 zmean   = 1;       
 zstdv   = .2;      
 [z,prob] = markovchain(ne,0.5,0.5,zstdv,zmean)

%% State-space S = ZxK
 k = linspace(0.001,0.5,500)';
 
 n = length(z)*length(k); 
 m = length(k);

 [Z,K] = gridmake(z,k);

%% Exact Solution (evaluated at zmean)
 kex = alpha*beta*zmean*K.^alpha    ;

%% Feasible Consumption set C>0 and Utility function
 C  = zeros(n,m);

 for i=1:m    
    C(:,i)=Z.*K.^alpha - k(i);   
 end
 if gamma==1;
    u = log(C);
 else
    u = (C.^(1-gamma))./(1-gamma); 
 end
 
 u(C<=0)=-Inf; 

%% Transition probability matrix 
 P = kron(speye(m,m),repmat(prob,m,1));

%% Bellman equation
 [v,x,pstar] = solvedp(u,P,beta,'policy');  clear P u C;

%% Stationary density
 d = markov(pstar);

c = reshape((Z.*K.^alpha - k(x)),500,ne);
kp = reshape(k(x),500,ne);
%% Plot figures 
 figure(4); bar(d)
 kmean = k(x)'*d
 cmean = zmean*kmean^alpha-kmean

% Numeric vs exact solution
 figure(5); scatter(K,k(x)); hold on; scatter(K,kex); hold off
