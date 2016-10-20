% MENDOZA_COLOMBIA  Mendoza (1991) Small Open Economy RBC Model 
%                   calibrated for the Colombian economy
%
% Written by F. Hamann. Feel free to copy, change and distribute
 close all; clear all; %clc;

 fprintf('\nMendoza 1991 model calibrated for Colombia\n')

%% Model Parameters (values taken from AER published paper)
 gamma  = 4.0;            % Risk aversion coefficient in PATACON  4.0
 omega  = 1.5;            % Inv. Frisch: 0.5 from Prada (2009)    1.5
 beta   = 0.095;          % to match 30% debt to GDP ratio        0.095
 r      = 0.04;           % foreign interest rate in PATACON     *0.035*
 delta  = 0.1;            % depreciation rate in PATACON          0.1
 alpha  = 0.31;           % capital's share of income in PATACON  0.31
 phi    = 0.024;          % Investment adj. cost to match I sdv   0.024
 rho    = 0.6;            % auto-correlation                      0.6
 rhoen  = 0;              % cross-correlation                     0
 sigmae = 0.0057*1.5;     % sdv TFP shocks to match GDP sdv       0.0057
 sigman = 0.0118;         % sdv R* shock taken from Mendoza       0.0118

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
 dsc = 1/(exp(-beta*log(1+css-((hss.^omega)./omega))))-1;

%% Construct state-space, reward function and transition matrix
% Number of actions and states 
 ne = 3; 
 n  = [ne*ne 33 33];                         % # of states; 
 m  = [n(2) n(3)];                           % # actions (endog. states)
 
 N = prod(n);                                % Total number of states
 M = prod(m);                                % Total number of actions

% Probability transition matrix
 disp('Discretizing VAR(1) as Markov Chain ... ')
 A0x    = [rho 0; 0 rho];
 vex    = [sigmae^2 0; 0 sigman^2];
 ntune  = 5000;
 type   = 1;
 [ps,s] = var_Markov_MM_General(A0x,vex,ne,ntune,type);
 P      = kron(speye(n(3)),kron(speye(n(2),n(2)),repmat(ps,n(2)*n(3),1)));
 disp('Done!'); 

% State space: SxAxK  
% Bounds for AxK

%  For beta = 0.095 and r = 0.035 -> ass<0
%   amin = ass*2.40 ; amax = -ass*0.40 ; a = linspace(amin,amax,n(2));
%   kmin = kss*0.96 ; kmax = kss*1.04  ; k = linspace(kmin,kmax,n(3));

%  For beta = 0.095 and r = 0.04 -> ass>0
% amin = ass*0.80 ; amax = ass*1.20 ; a = linspace(amin,amax,n(2));
% kmin = kss*0.97 ; kmax = kss*1.03 ; k = linspace(kmin,kmax,n(3));
% For wider [A K] bounds, set n = [4 44 44] to match GDP sdv.
  amin = ass*0.70 ; amax = ass*1.30 ; a = linspace(amin,amax,n(2));
  kmin = kss*0.95 ; kmax = kss*1.05 ; k = linspace(kmin,kmax,n(3));
 
 [Se,Sn,A,K] = gridmake(exp(s),a',k');
 [aa,kk]     = gridmake(a',k');

% Construct the reward function, f
 c = zeros(N,M);
 L = ((1-alpha)*Se.*K.^alpha).^(1/(alpha+omega-1));
 Y = Se.*K.^alpha.*L.^(1-alpha);

 for i=1:M
  c(:,i)=(Y+(1-delta)*K+(1+r*Sn).*A-aa(i)-kk(i)-(phi/2)*(kk(i)-K).^2); 
 end

 L = repmat(L,1,M); % replicate matrix to be conformable with c

 if gamma == 1
  u = log(c)-3*L; 
 else
  u = ((c-((L.^omega)./omega)).^(1-gamma)-1)/(1-gamma);
 end
 u(c<=0) = NaN;                                 

%% Solve the Model

 betahat = exp(-beta*log(1+c-((L.^omega)./omega)));
 betahat(betahat>=1) = NaN;

 disp('Solving the model iterating on the Bellman equation...'); 
 [v,x,pstar] = solvedp(u,P,betahat);  
 disp('Done!'); 

 ax = aa(x)';                          % Store a' policy in ax vector
 kx = kk(x)';                          % Store k' policy in ax vector
 lx = ((1-alpha)*kk(x).^alpha).^(1/(alpha+omega-1));
 cx = (Y+(1-delta)*K+(1+r*Sn).*A-aa(x)-kk(x)-(phi/2)*(kk(x)-K).^2); 

%% Ergodic moments 

 pi    = ergdist(pstar);                           
 kmean = pi'*kk(x);                                  
 amean = pi'*aa(x);                                 
 lmean = pi'*((1-alpha)*kk(x).^alpha).^(1/(alpha+omega-1));
 ymean = pi'*(kk(x).^alpha.*lx.^(1-alpha));
 cmean = pi'*cx;
 cmean = ymean+r*amean-delta*kmean;
 imean = delta*kmean;
 d2yss = amean/ymean;

%% Simulated moments

 T      = 10000;      
 s0     = getindex([amean kmean],[A K]);   
 spath  = simulmarkov(pstar,T,s0);

 Lpath  = ((1-alpha)*Se(spath).*K(spath).^alpha).^(1/(alpha+omega-1));
 ypath  = Se(spath).*K(spath).^alpha.*Lpath.^(1-alpha);
 ipath  = kx(spath)' - (1-delta)*K(spath);
 CApath = aa(x(spath))-A(spath);
 ACpath = (phi/2)*(kk(x(spath))-K(spath)).^2;
 cpath  = ypath - ipath + r.*Sn(spath).*A(spath)-CApath-ACpath;
 
%% Graphics 
% Page   : indexes to states of productivity (last index): 1=high, 2=low
% Row    : indexes to states of capital (first index)
% Column : indexes to states of asset (second index)

 PP = reshape(pi,n(1),n(2),n(3));          % Rearrange pi 
 PP = permute(PP,[3 2 1]);                 % Flip row and page subscripts

 figure('Color',[1 1 1]);                  % White background color
 xlabel('K'),ylabel('A'),zlabel('PKA')
 surf(k,a,PP(:,:,1)','FaceColor',[1 1 1])

%% Simulation
 tau = 3;
 [mn,sd,corr,acorr] = samplemoms([log(ypath) log(cpath) log(ipath) ...
                                  log(Lpath) log(K(spath)) log(A(spath)) ...
                                  CApath./ypath],1,tau);

%% Reporting main moments

fprintf('\n\n\n')
fprintf('                     Small Open Economy RBC Model         \n')
fprintf('                 calibrated for the Colombian economy     \n')
fprintf('                      Great Ratios (pct of GDP)           \n')
fprintf('\nRatio             Colombia  Det.SS   Ergodic Simulated  Level')
fprintf('\nOutput           %8.2f %8.2f %8.2f %8.2f %8.3f'   ,1,yss/yss,ymean/ymean,exp(mn(1))/exp(mn(1)),exp(mn(1)))
fprintf('\nConsumption      %8.2f %8.2f %8.2f %8.2f %8.3f'   ,0.786,css/yss,cmean/ymean,exp(mn(2))/exp(mn(1)),exp(mn(2)))
fprintf('\nInvestment       %8.2f %8.2f %8.2f %8.2f %8.3f'   ,0.254,iss/yss,imean/ymean,exp(mn(3))/exp(mn(1)),exp(mn(3)))
fprintf('\nLabor            %8.2f %8.2f %8.2f %8.2f %8.3f'   ,NaN,hss,lmean,exp(mn(4)),lmean)
fprintf('\nK/Y              %8.2f %8.2f %8.2f %8.2f %8.3f'   ,2,kss/yss,kmean/ymean,exp(mn(5))/exp(mn(1)),exp(mn(5)))
fprintf('\nNFA to output    %8.2f %8.2f %8.2f %8.2f %8.3f\n' ,-0.3,ass/yss,amean/ymean,exp(mn(6))/exp(mn(1)),exp(mn(6)))

fprintf('\n')
fprintf('                                  Moments                    \n')
fprintf('                            Data                        Model  ')
fprintf('\nVariable          Std.Dev  Corr(i,y) Autocorr  Std.Dev Corr(i,y) Autocorr')
fprintf('\nOutput           %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,1.58,1.00,0.76,sd(1),corr(1,tau+1),acorr(tau,1))
fprintf('\nConsumption      %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,1.50,0.83,0.80,sd(2),corr(2,tau+1),acorr(tau,2))
fprintf('\nInvestment       %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,8.34,0.65,0.81,sd(3),corr(3,tau+1),acorr(tau,3))
fprintf('\nNFA              %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,NaN ,NaN ,NaN ,sd(6),corr(6,tau+1),acorr(tau,6))
fprintf('\nCurrent account  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' ,NaN ,NaN ,NaN ,sd(7),corr(7,tau+1),acorr(tau,7))

