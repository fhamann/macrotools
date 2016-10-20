% BROCK_MIRMAN Brock and Mirman (1972) optimal stochastic growth model 
%
% Written by F. Hamann. Feel free to copy, change and distribute
 fprintf('\nBrock and Mirman model \n')

%% Parameters
 gamma	= 1;        % risk aversion 
 beta	= 0.98;     % discount factor
 alpha  = 0.3;      % capital share

%% Markov chain for z
 zmean   = 1;       
 zstdv   = .1;      
 [z,prob] = markovchain(5,0.5,0.5,zstdv,zmean);

%% State-space S = ZxK
 k = linspace(0.001,0.5,500)';
 
 n = length(z)*length(k); 
 m = length(k);

 [Z,K] = gridmake(z,k);

%% Exact Solution (evaluated at zmean)
 kex = alpha*beta*zmean*K.^alpha;    

%% Conjunto de consumo factible C>0 y funcion de utilidad
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

%% Plot figures 
 figure(1); bar(d)
 kmean = k(x)'*d
 cmean = zmean*kmean^alpha-kmean

% Numeric vs exact solution
 figure(2); scatter(K,k(x)); hold on; scatter(K,kex); hold off