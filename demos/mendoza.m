% MENDOZA Mendoza (1991) Small Open Economy Real Business Cycle Model 
%
% Written by F. Hamann. Feel free to copy, change and distribute
 fprintf('\nMendoza 1991 model \n')

%% Model Parameters (values taken from AER published paper)
 gamma  = 2.0;                             % Risk aversion coefficient 
 omega  = 1.455;                           % Labor elasticity
 beta   = 0.11;                            % discount factor constant  
 r      = 0.04;                            % interest rate 
 delta  = 0.1;                             % depreciation
 alpha  = 0.32;                            % capital's share of income
 phi    = 0.028;                           % adjustment costs 0.025
 rho    = 0.42;                            % auto-correlation rho=0.356;
 rhoen  = rho;                             % cross-correlation 
 sigmae = 0.0118;                          % stdev of productivity shocks
 sigman = 0.0118;                          % stdev of interest rate shock

%% Deterministic steady state (sometimes is useful)
z   = 0; 
R   = (1+r);
hss = ((1-alpha) * (alpha/(r+delta))^(alpha/(1-alpha)))^(1/(omega-1));
css = R^(1/beta) + hss^omega/omega -1; 
kss = hss*((r+delta)/alpha)^(-1/(1-alpha)); 
iss = delta*kss;
yss = exp(z)*kss^alpha*hss^(1-alpha);
tb  = yss - css - iss;
ass = -tb/r;
tby = tb/yss;
cay = (r*ass + tb)/yss; 
a2y = -ass/yss;

%% Construct state-space, reward function and transition matrix
% Number of actions and states 
 n = [4 22 22];                              % # states n = [4 22 22]; 
 m = [n(2) n(3)];                            % # actions (endog. states)
 
 N = prod(n);                                % Total number of states
 M = prod(m);                                % Total number of actions

% Probability transition matrix
 prob = ptm(rho,rhoen);
 PI   = repmat(prob,n(2)*n(3),1);                         
 P    = kron(speye(n(3)),kron(speye(n(2),n(2)),PI));

% State space: SxAxK 
 s1   = exp([-sigmae;sigmae]); 
 s2   = exp([-sigman;sigman]);
 
% Original bounds in Mendoza(1991) 
%  kmin = 3.25;  kmax = 3.56;  k = linspace(kmin,kmax,n(3));               
%  amin = -1.42; amax = 0.08;  a = linspace(amin,amax,n(2));               

 kmin = kss*0.9565; kmax = kss*1.0478; k = linspace(kmin,kmax,n(3));
 amin = ass*2.4060; amax = ass*0.1356; a = linspace(amin,amax,n(2));

 [S1,S2,A,K] = gridmake(s1,s2,a',k');           
 [aa,kk]     = gridmake(a',k');

% Construct the reward function, f
 c = zeros(N,M);
 L = ((1-alpha)*S1.*K.^alpha).^(1/(alpha+omega-1));
 Y = S1.*K.^alpha.*L.^(1-alpha);

for i=1:M
  c(:,i)=(Y+(1-delta)*K+(1+r*S2).*A-aa(i)-kk(i)-(phi/2)*(kk(i)-K).^2); 
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
 ysdv  = sqrt(pi'*(S1.*kk(x).^alpha.*lx.^(1-alpha)-ymean).^2);

 cmean = ymean+r*amean-delta*kmean;
 c     = (Y+(1-delta)*K+(1+r*S2).*A-aa(x)-kk(x)-(phi/2)*(kk(x)-K).^2); 
 csdv  = sqrt(pi'*(c-cmean).^2);

 d2yss = amean/ymean;

%% Graphics 
 ax = aa(x)';                                % Store policy in ax vector
 kx = kk(x)';
 V  = reshape(v,n(1),n(2),n(3));             % Rearrange v 
 X  = reshape(x,n(1),n(2),n(3));             % Rearrange x
 PP = reshape(pi,n(1),n(2),n(3));            % Rearrange pi 
 Ax = reshape(ax,n(1),n(2),n(3));            % Rearrange a(x) in a matrix
 Kx = reshape(kx,n(1),n(2),n(3));            % Rearrange k(x) in a matrix

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

 figure3=figure('Color',[1 1 1]); % color de fondo: blanco
% xlim([kmin kmax]);
% ylim([amin amax]);
% zlim([0 0.015]);
 xlabel('K'),ylabel('A'),zlabel('PKA')
 surf(k,a,PP(:,:,1)','FaceColor',[1 1 1])


%% Statistics
% Reporting main theoretical moments

names   = ['Output               ';
           'Consumption          ';    
           'Capital              ';
           'Labor                ';
           'Net foreign assets   ';
           'NFA to GDP           '];

ss      = [yss;   css;   kss;   hss;   ass;   ass/yss];       
means   = [ymean; cmean; kmean; lmean; amean; amean/ymean];
espacio = ['   ';'   ';'   ';'   ';'   ';'   '];
stddevs = [ysdv;csdv;ksdv;lsdv;asdv; NaN];

display('Steady state, stoch. mean and std dev. of Mendoza (1991)')
TABLE=[names num2str(ss) espacio num2str(means) espacio num2str(stddevs)]