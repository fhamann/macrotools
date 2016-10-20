% RISK_MATTERS  Risk Matters: The Real Effects of Volatility Shocks
%               Jesús Fernández-Villaverde, Pablo Guerrón-Quintana, 
%               Juan F. Rubio-Ramírez, and Martin Uribe
%               American Economic Review 101 (October 2011): 2530-2561

% Written by F. Hamann. Feel free to copy, change and distribute

 mendoza; vini = v;                         % runs Mendoza (1991) SOE
beta   = 1/(1+r);                          % discount factor constant  
phik   = phi;                              % adjustment costs 0.025
phia   = 0.00074;                          % no adjustment costs, phi=0
Abar   = amean;                            % Steady state in Mendoza


% Model Parameters (taken from Uribe & Schmitt-Grohe
%  gamma  = 2.0;                              % Risk aversion coefficient 
%  omega  = 1.455;                            % Labor elasticity
%  delta  = 0.1;                              % depreciation
%  alpha  = 0.32;                             % capital share
%  rho    = 0.42;                             % shocks auto-correlation 
%  rhoen  = 0.42;                             % shocks cross-correlation 
%  sigmae = 0.0129;                           % stdev of productivity shock
%  sigman = 0.0129;                           % stdev of interest rate shoc
%  beta   = 0.96;                             % discount factor constant  
%  phik   = 0.0280;                           % adjustment costs 0.025
%  phia   = 0.00074;                          % no adjustment costs, phi=0
%  Abar   = -0.7442;                          % Steady state in Mendoza

% Number of actions and states 
%  n = [4 22 22];                             % # states n = [4 22 22]; 
%  m = [n(2) n(3)];                           % # actions (endog. states)
% 
%  N = prod(n);                               % Total number of states
%  M = prod(m);                               % Total number of actions
% 
% %% Construct state-space, reward function and transition matrix
% % Probability transition matrix
%  prob = ptm(rho,rhoen);
%  PI   = repmat(prob,n(2)*n(3),1);                         
%  P    = kron(speye(n(3)),kron(speye(n(2),n(2)),PI));
% 
% % State space: SxAxK 
%  s1   = exp([-sigmae;sigmae]);
%  s2   = exp([-sigman;sigman]);
% 
%  kmin = 3.25;  kmax = 3.56;  k = linspace(kmin,kmax,n(3));               
%  amin = -1.42; amax = 0.08;  a = linspace(amin,amax,n(2));               
% 
%  [S1,S2,A,K] = gridmake(s1,s2,a',k');           
%  [aa,kk] = gridmake(a',k');

% Reward function, u 
  c = zeros(N,M);
  f = zeros(N,M);
  L = ((1-alpha)*S1.*K.^alpha).^(1/(alpha+omega-1));
  Y = S1.*K.^alpha.*L.^(1-alpha);

 for i=1:M
   c(:,i) = (Y+(1-delta)*K+(1+r*S2).*A-aa(i)-kk(i)...
            -(phik/2)*(kk(i)-K).^2)-(phia/2)*(aa(i)-Abar).^2; 
 end

 L = repmat(L,1,M); % replicate matrix to be conformable with c

 if gamma == 1
   u = log(c)-3*L; 
 else
   u = ((c-((L.^omega)./omega)).^(1-gamma)-1)/(1-gamma);
 end
 
 u(c<=0) = NaN;

%% Solve the Model
 [v,x,pstar] = solvedp(u,P,beta,'value',vini);   

%% Ergodic moments (means and std deviations)

 pi = ergdist(pstar);                           

 kmean = pi'*kk(x);                                  
 ksdv  = sqrt(pi'*(kk(x)-kmean).^2);                       

 amean = pi'*aa(x);                                 
 asdv  = sqrt(pi'*(aa(x)-amean).^2);                       

 lmean = ((1-alpha)*kmean.^alpha).^(1/(alpha+omega-1));
 lx    = ((1-alpha)*kk(x).^alpha).^(1/(alpha+omega-1));
 lsdv  = sqrt(pi'*(lx-lmean).^2);

 ymean = kmean^alpha*lmean^(1-alpha);
 ysdv  = sqrt(pi'*(S1.*kk(x).^alpha.*lx.^(1-alpha)-ymean).^2);

 cmean = ymean+r*amean-delta*kmean;
 c     = (Y+(1-delta)*K+(1+r*S2).*A-aa(x)-kk(x)...
         -(phi/2)*(kk(x)-K).^2-(phia/2)*(aa(x)-Abar).^2); 
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

 figure4=figure('Color',[1 1 1]); % color de fondo, en este caso es blanco
 axes1 = axes(...
  'CameraPosition',[4.913 -9.243 0.05396],...
  'CameraUpVector',[-15.78 89.04 0.7743],...
  'Parent',figure4);             % configurando el eje
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

meansfv = [ymean;cmean;kmean;lmean;amean;amean/ymean];
espacio = ['   ';'   ';'   ';'   ';'   ';'   '];
stddevsfv = [ysdv;csdv;ksdv;lsdv;asdv; NaN];

display('Means & Std.Devs of Mendoza vs Fernandez-Villaverde et al.')
TABLE = [names num2str(means)   espacio num2str(meansfv)  espacio...
               num2str(stddevs) espacio num2str(stddevsfv) ]
