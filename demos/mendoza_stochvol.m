% MENDOZA_STOCHVOL Mendoza (1991) Small Open Economy RBC Model with
%                  stochastic volatility in interest rate process

% Written by F. Hamann. Feel free to copy, change and distribute

 clear all

% Model Parameters
 gamma  = 2.0;                               % Risk aversion coefficient 
 omega  = 1.455;                             % Labor elasticity
 beta   = 0.11;                              % discount factor constant  
 r      = 0.04;                              % interest rate 
 delta  = 0.1;                               % depreciation
 alpha  = 0.32;                              % capital's share of income
 phi    = 0.028;                             % adjustment costs 0.025
 
%% Construct state-space, reward function and transition matrix
% Number of actions and states 
 n = [4 44 44];                              % # states n = [4 22 22]; 
 m = [n(2) n(3)];                            % # actions (endog. states)
 
 N = prod(n);                                % Total number of states
 M = prod(m);                                % Total number of actions

% Probability transition matrix (stochastic volatility Markov process)

 p      = 0.8;         % probability of staying in low interest rates
 q      = 0.8;         % probability of staying in high interest rates
 pL     = 0.9;         % probability of staying in low volatility regime
 pH     = 0.8;         % probability of staying in high volatility regime
 Rm     = 1;           % interest rate mean
 sigman = 0.0118;      % stdev of interest rate shock
 mean   = 0;           % mean of shocks (in logs, later exp -> 1)
 exvol  = 1;           % excess volatility in tails
 
 [s,prob] = svmc(p,q,pL,pH,sigman,mean,exvol);  s = exp(s);
  
 PI   = repmat(prob,n(2)*n(3),1);                         
 P    = kron(speye(n(3)),kron(speye(n(2),n(2)),PI));

% State space: SxAxK 

 kmin = 3.25;  kmax = 3.56;  k = linspace(kmin,kmax,n(3));              
% kmin = 3.2;  kmax = 3.6;  k = linspace(kmin,kmax,n(3));  % for exvol>=10            
 amin = -1.42; amax = 0.08;  a = linspace(amin,amax,n(2));               

 [S,A,K] = gridmake(s,a',k');           
 [aa,kk] = gridmake(a',k');

% Construct the reward function, u
 c = zeros(N,M);
 L = ((1-alpha)*K.^alpha).^(1/(alpha+omega-1));
 Y = K.^alpha.*L.^(1-alpha);

for i=1:M
  c(:,i)=(Y+(1-delta)*K+(1+r*S).*A-aa(i)-kk(i)-(phi/2)*(kk(i)-K).^2); 
end

L = repmat(L,1,M); % replicate matrix to be conformable with c

if gamma == 1
  u = log(c)-3*L; 
else
  u = ((c-((L.^omega)./omega)).^(1-gamma)-1)/(1-gamma);
end
 u(c<=0) = NaN;                                 

% Construct the transition matrix
 PI = repmat(prob,n(2)*n(3),1);                         
 I  = speye(n(2),n(2));                         
 P  = kron(speye(n(3)),kron(I,PI));

% Solve the Model
 betahat = exp(-beta*log(1+c-((L.^omega)./omega)));
 betahat(betahat>=1) = NaN;

 [v,x,pstar] = solvedp(u,P,betahat);  % Value, Policy and Transition

 r*s
%% Ergodic moments 

 pi = ergdist(pstar);                           

 kmean = pi'*kk(x);                                  
 ksdv  = sqrt(pi'*(kk(x)-kmean).^2);                       

 amean = pi'*aa(x);                                 
 asdv  = sqrt(pi'*(aa(x)-amean).^2);                       

 lmean = ((1-alpha)*kmean.^alpha).^(1/(alpha+omega-1));
 lx    = ((1-alpha)*kk(x).^alpha).^(1/(alpha+omega-1));
 lsdv  = sqrt(pi'*(lx-lmean).^2);

 ymean = kmean^alpha*lmean^(1-alpha);
 ysdv  = sqrt(pi'*(kk(x).^alpha.*lx.^(1-alpha)-ymean).^2);

 cmean = ymean+r*amean-delta*kmean;
 c     = (Y+(1-delta)*K+(1+r*S).*A-aa(x)-kk(x)-(phi/2)*(kk(x)-K).^2); 
 csdv  = sqrt(pi'*(c-cmean).^2);

 d2yss = amean/ymean;

%% Graphics 
 ax = aa(x)';                                  % Store policy in ax vector
 kx = kk(x)';
 V  = reshape(v,n(1),n(2),n(3));               % Rearrange v 
 X  = reshape(x,n(1),n(2),n(3));               % Rearrange x
 PP = reshape(pi,n(1),n(2),n(3));              % Rearrange pi 
 Ax = reshape(ax,n(1),n(2),n(3));              % Rearrange a(x) in a matrix
 Kx = reshape(kx,n(1),n(2),n(3));              % Rearrange k(x) in a matrix

 V = permute(V,[3 2 1]);    % Reverse row and page subscripts 
 X = permute(X,[3 2 1]);    % Reverse row and page subscripts
 PP = permute(PP,[3 2 1]);  % Reverse row and page subscripts
 Ax = permute(Ax,[3 2 1]);  % Reverse row and page subscripts
 Kx = permute(Kx,[3 2 1]);  % Reverse row and page subscripts

% Page   : indexes to states of productivity (last index): 1=high, 2=low
% Row    : indexes to states of capital (first index)
% Column : indexes to states of asset (second index)

% figure(1)
% surf(k,a,V(:,:,1)')
% 
% figure(2)
% plot(aa(x));

 figure3=figure('Color',[1 1 1]); % color de fondo, en este caso es blanco
 axes1 = axes(...
  'CameraPosition',[4.913 -9.243 0.05396],...
  'CameraUpVector',[-15.78 89.04 0.7743],...
  'Parent',figure3);             
 grid(axes1,'off')
 xlim(axes1,[kmin kmax]);
 ylim(axes1,[amin amax]);
 zlim(axes1,[0 0.015]);
 xlabel('K'),ylabel(axes1,'A'),zlabel(axes1,'PKA'),hold(axes1)
 surf(k,a,PP(:,:,1)','FaceColor',[1 1 1])
%toc

%% Statistics
% Reporting main theoretical moments

names   = ['Output               ';
           'Consumption          ';    
           'Capital              ';
           'Labor                ';
           'Net foreign assets   ';
           'NFA to GDP           '];

means = [ymean;cmean;kmean;lmean;amean;amean/ymean];
espacio = ['   ';'   ';'   ';'   ';'   ';'   '];
stddevs = [ysdv;csdv;ksdv;lsdv;asdv; NaN];

display('Means and Std Deviations of Mendoza (1991) ')
TABLE = [names num2str(means) espacio num2str(stddevs)]
