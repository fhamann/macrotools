%% BAH_UE.M  Bewley-Aiyagari-Huggett vs. Uzawa-Epstein model
%            Follows Durdu, Mendoza and Terrones (2007) NBER-WP 13123
%
% Written by F. Hamann. Feel free to copy, change and distribute.

% % % % Features
% % 1) Parameters:
% It is included:the share of imported inputs in gross product
%                the VAR(1) coefficient matrix for the discretization as Markov chain             
% % 2) Normalization and data to match:
% It is included: gpN - Non-Tradables Gross-Product (GDP Non-Tradables plus imported inputs)
%                 imp - Imported inputs  
%                  Z  - TPF trend mean
% It is changed:  AN  - NT Absorption: NT_Gross_Product - NT_consumption 
% % 3) Markov Chain approximation:
% It is included: Discretization of a VAR(1) as a Markov chain
%                    + Gospondinov, N. & Lkhagvasuren, D. (2014) approach.
%                      - article: DOI: 10.1002/jae.2354
%                      - code: https://sites.google.com/site/dlkhagva/var_mmm.
%                    + Tauchen(1986) approach. Coded by Iskander Karibzhanov
%                      - article https://doi.org/10.1016/0165-1765(86)90168-0
%                      - code: available in http://karibzhanov.com/
% % 4) State space for b:
% The shocks for yT and gpN are included and it is created a finer grid
% % 5) Utility function and feasible consumption
% It is included now an equation for the consumption of NT using the new grid
% and the imported imputs in the Tradables consumption as a function of GPN
% (alpha*GPN)


% % % Note
% It is worth mentioning that using both, Gospondinov, N. &
% Lkhagvasuren, D. (2014) and Tauchen(1986) approaches, the VAR(1)
% discretization as a Markov chain does not match the unconditional moments
% of the Mexican data for period 1965-2005 using the provided coefficient
% matrix and the error terms variance-covariance matrix in the article.



 fprintf('\nDurdu, Mendoza and Terrones (2007) NBER-WP 13123')
 fprintf('\nSingle good SOE model calibrated for Mexico \n')

%% Parameters
 sigma = 2;                         % risk aversion coefficient 2
 R     = 1.059;                     % gross asset return rate1.059
 mu    = 0.316;                     % elasticity of substituion
 zss   = 0.98;                      % productivity shock  
 alpha = 0.2;                       % Share of imported inputs in gross product
 rhoyT = 1.088;                     % First order lag in yT VAR(1) equation respect to yT .108
 rhoyN = 0.3;                       % First order lag in yN VAR(1) equation respect to yN .3
 rhoyTN= 0.564;                     % First order lag in yT VAR(1) equation respect to yN .564
 rhoyNT= -0.655;                    % First order lag in yT VAR(1) equation respect to yT .655


%% Normalization and data to match
% Set pN=1 & yT+pN*yN = 1 to get SS alloc. as % of GDP
 pN      = 1;
 pm      = 1;
 pNyN2yT = 1.543;                    % NT to T output ratio 1.543 
 bss     = -0.44;                    % NFA to GDP ratio     -0.44

 phi    = 0.71;                     % ad-hoc borrowing limit    -0.71
 
 yT     = 1/(1+pNyN2yT);             % T output
 yN     = 1-yT;                      % NT value of output (in T-units)


 cT2yT     = 0.665;                  % Tradables Consumption to Output ratio
 pNcN2yN   = 0.71;                   % Non-Tradables Consumption to Output ratio
 cT        = cT2yT*yT;               % Tradables Consumption
 cN        = pNcN2yN*(yN/pN);        % Non Tradables Consumption
 a         = 1/(1+(pN/(cT/cN)^(1+mu))); % CES weighting factor
 gpN       = yN/(1-alpha);           % Non-Tradables Gross-Product nontradables (GDP plus imported inputs)
 imp       = alpha*gpN ;             % Imported inputs
 AT        = yT+bss*(R-1)-cT-imp;    % Tradables Absorption
 AN        = gpN-cN;                 % Non-Tradables Absorption
 Z         = (gpN^(1-alpha))/(alpha^alpha); % TPF trend mean


%% Discount factor
 beta_bah = 0.94395;                    % Bewley-Aiyagari-Huggett   0.94
 css = (a*(cT^(-mu))+(1-a)*(cN^(-mu)))^(-1/mu);
 disc_ue  = log(R)/log(1+css);        % Uzawa-Epstein disc rate   
 beta_ue  = (1+css)^-disc_ue;         % Uzawa-Epstein disc factor

 
%% Stationarity condition
 if R*beta_bah>=1; display('Set beta*R<1 for convergence'); end;

%% Markov chain for y
%  nz   = 3;
%  rhoz = 0.597;                       % GDP cycle autocorr.  0.597
%  sdvz = 0.03301;                     % GDP cycle std dev.   0.03301
%  sdve = sqrt(sdvz^2*(1-rhoz^2))
%  p    = (1+rhoz)/2;
%  q    = p;
% 
%  [z,Pz] = markovchain(nz,p,q,sdvz,zss);
%  [z,Pz] = markovchain(nz,p,q,sdvz,Z);
 
%  spath  = simulmarkov(Pz,50000,1); zpath  = z(spath);
%  
%  [zmean,zsdv,zcorrcont] =samplemoms(zpath,1,2);




%% VAR(1) Markov Chain approximation
 disp('Discretizing VAR(1) as Markov Chain ... ')
 ne     = 3;
A0x    = [rhoyT rhoyTN; rhoyNT rhoyN];      % Coefficient Matrix


% % % Gosponidov & Lhagkvasuren approach 
%  vex    = [0.000601 0.00055; 0.00055 0.0012];  
%  ntune  = 10000;
%  type   = 0;
%  [Pz,z] = var_Markov_MM_General(A0x,vex,ne,ntune,type);

 % % Tauchen Approach
nn = [3 3];
mumc = [yT yN]; 
covae = [0.000601 0.0012];
vee    = [sqrt(covae(1,1)) sqrt(covae(1,2))]; 
[z,Pz] = tauchen(nn,mumc,A0x,vee);
 
 
% % % % 
z = z';
Pz = Pz';
T = 2000000;  % T
T_0 = 200000; % Not-usable Data 
S = simulmarkov(Pz,T_0+T);
S_efect=S(T_0+1:end);
zsim=z(S_efect,:);

zmean=mean(zsim)
varzsim=cov(zsim)
stdzsim = diag([sqrt(varzsim(1,1)) sqrt(varzsim(2,2))]) 
rhoze=varzsim(1,2)/sqrt(varzsim(1,1)*varzsim(2,2))
autcorrT=autocorr(zsim(:,1),1)
autcorrN=autocorr(zsim(:,2),1)


%% State-space S = YxB
b = linspace(-phi,phi,450)';   

 [YT,GPN,B] = gridmake(z(:,1)*yT,z(:,2)*gpN,b);  % For the VAR(1) approach
 
 %  [YT,GPN,B] = gridmake(z*yT,z*gpN,b);         % For the AR(1) approach

 n = length(z)*length(z)*length(b);
 m = length(b);

%% Utility function and feasible consumption C>=0
 CT = ones(n,m);
 CN = ones(n,m); 
 
 for i=1:m    
    CN(:,i)= GPN+AN;    
    CT(:,i)= YT+R*B+AT-b(i)-alpha*GPN;
    end

   CT(CT<=0) = NaN;
   CN(CN<=0) = NaN;
   c    = (a*(CT.^(-mu)) + (1-a)*(CN.^(-mu))).^(-1/mu); 
   u    = (c.^(1-sigma))./(1-sigma); 
   u(c<=0)=-Inf;
 
 beta = exp(-disc_ue*log(1+c));
 beta(beta>=1) = NaN;

%% Transition probability matrix
 P = kron(speye(ne*ne,ne*ne),kron(speye(m,m),repmat(Pz,m,1))); % VAR(1)

%  P = kron(speye(nz,nz),kron(speye(m,m),repmat(Pz,m,1)));       % AR(1)

%% Solve Bellman equation
 [v,x,pstar] = solvedp(u,P,beta_bah,'policy');  % Bewley-Aiyagari-Huggett

% [v,x,pstar] = solvedp(u,P,beta);                % Uzawa-Epstein
 clear P u C;

%% Steady State Distribution
 d = ergdist(pstar);

%% Summary statistics 
  cTx = YT+R*B-b(x)+AT-alpha*GPN; 
  cNx = GPN+AN;
  pN = (1-a)/a*(cTx./cNx).^(1+mu);

 zmean1 = ergdist(Pz)'*z(:,1);
 zmean2 = ergdist(Pz)'*z(:,2);
 bmean = b(x)'*d;
 cmean = c'*d;

%% Plot some model properties
%  plotdp(v,x,pstar,YT,B,z*yT,b);



%% Simulation
 T      = 500;      
 s0     = findnearest(bmean,B);   
 spath  = simulmarkov(pstar,T,s0);

 yTpath  = YT(spath);
 cTpath  = yTpath + R*B(spath)-b(x(spath))+ AT(spath)-alpha*GPN(spath);
 cNpath  = GPN(spath) + AN(spath);
 pNpath  = (1-a)/a*(cTpath./cNpath).^(1+mu);
 CApath  = b(x(spath))-B(spath);

 figure(2)
 plot([yTpath cTpath])

 tau =3;
 [mn,sd,corr,acorr,auocorr] = samplemoms([yTpath cTpath B(spath) CApath 1./pNpath cNpath],1,tau);

%% Model steady state statistics  
fprintf('\n\n\n')
fprintf('                     Small Open Economy RBC Model         \n')
fprintf('                  calibrated for the Mexican economy     \n')
fprintf('                      Great Ratios (pct of GDP)           \n')
fprintf('\nRatio             Mexico  Det.SS   Ergodic Simulated  Level')
fprintf('\nOutput           %8.2f %8.2f %8.2f %8.2f %8.2f'   ,1,yT,yT/yT,mn(1)/mn(1),mn(1))
fprintf('\nConsumption T    %8.2f %8.2f %8.2f %8.2f %8.2f'   ,0.786,cT/yT,cT/yT,mn(2)/mn(1),mn(2))
fprintf('\nNFA to output    %8.2f %8.2f %8.2f %8.2f %8.2f'   ,bss  ,bss/yT,bmean/yT,mn(3)/mn(1),mn(3))
fprintf('\nReal Exch. rate  %8.2f %8.2f %8.2f %8.2f %8.2f'   ,bss  ,bss/yT,bmean/yT,mn(5),mn(5))
fprintf('\nTime b.c.        %8.2f %8.2f %8.2f %8.2f %8.2f\n' ,0.162,0,sum(d(1:nz)),NaN,NaN)

fprintf('\n')
fprintf('                                  Moments                    \n')
fprintf('                            Data                        Model  ')
fprintf('\nVariable          Std.Dev  Corr(i,y) Autocorr  Std.Dev Corr(i,y) Autocorr')
fprintf('\nOutput           %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,2.47,1.00,0.76,sd(1),corr(1,tau+1),acorr(tau,1))
fprintf('\nConsumption      %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,3.31,0.83,0.80,sd(2),corr(2,tau+1),acorr(tau,2))
fprintf('\nNFA              %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,NaN ,NaN ,NaN ,sd(3),corr(3,tau+1),acorr(tau,3))
fprintf('\nCurrent account  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,NaN ,NaN ,NaN ,sd(4),corr(4,tau+1),acorr(tau,4))
fprintf('\nReal Exch. rate  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' ,NaN ,NaN ,NaN ,sd(5),corr(5,tau+1),acorr(tau,5))




