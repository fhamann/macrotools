%% BAH_UE.M  Bewley-Aiyagari-Huggett vs. Uzawa-Epstein model
%            Follows Durdu, Mendoza and Terrones (2007) NBER-WP 13123
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nDurdu, Mendoza and Terrones (2007) NBER-WP 13123')
 fprintf('\nTwo-good SOE model calibrated for Mexico \n')

%% Parameters
 sigma = 2;                         % risk aversion coefficient 2
 R     = 1.059;                     % gross asset return rate1.059
 mu    = 0.316;                     % elasticity of substituion
 zss   = 0.98;                      % productivity shock  
 
%% Normalization and data to match
% Set pN=1 & yT+pN*yN = 1 to get SS alloc. as % of GDP
 pN      = 1;
 pNyN2yT = 1.543;                    % NT to T output ratio 1.543 
 bss     = -0.44;                    % NFA to GDP ratio     -0.44

 phi    = -0.71;                     % ad-hoc borrowing limit    -0.51
 
 yT     = 1/(1+pNyN2yT);             % T output
 yN     = 1-yT;                      % NT value of output (in T-units)

 cT2yT     = 0.665;
 pNcN2pNyN = 0.71;
 cT        = cT2yT*yT;
 cN        = pNcN2pNyN*(pN*yN);
 a         = 1/(1+(pN/(cT/cN)^(1+mu)));
 AT        = yT+bss*(R-1)-cT;
 AN        = yN-cN;

%% Discount factor
 beta_bah = 0.94395;                    % Bewley-Aiyagari-Huggett   0.94
 css = (a*(cT^(-mu))+(1-a)*(cN^(-mu)))^(-1/mu);
 disc_ue  = log(R)/log(1+css);        % Uzawa-Epstein disc rate   
 beta_ue  = (1+css)^-disc_ue;         % Uzawa-Epstein disc factor

%% Stationarity condition
 if R*beta_bah>=1; display('Set beta*R<1 for convergence'); break; end;

%% Markov chain for y
 nz   = 3;
 rhoz = 0.597;                       % GDP cycle autocorr.  0.597
 sdvz = 0.03301;                     % GDP cycle std dev.   0.03301
 sdve = sqrt(sdvz^2*(1-rhoz^2))
 p    = (1+rhoz)/2;
 q    = p;

 [z,Pz] = markovchain(nz,p,q,sdvz,zss);

 spath  = simulmarkov(Pz,50000,1); zpath  = z(spath);
 
 [zmean,zsdv,zcorrcont] =samplemoms(zpath,1,2);
 
%% State-space S = YxB
 b = linspace(phi,0,1000)';   
 
 [YT,B] = gridmake(z*yT,b);

 n = length(z)*length(b); 
 m = length(b);

%% Utility function and feasible consumption C>=0
 CT = cT*ones(n,m);
 CN = cN*ones(n,m); 
 
 for i=1:m    
    CT(:,i)=YT+R*B-b(i)-AT;  
 end

   CT(CT<=0) = NaN;
   c    = (a*(CT.^(-mu)) + (1-a)*(CN.^(-mu))).^(-1/mu); 
   u    = (c.^(1-sigma)-1)./(1-sigma); u(c<=0)=-Inf;
 
 beta = exp(-disc_ue*log(1+c));
 beta(beta>=1) = NaN;

%% Transition probability matrix
 P = kron(speye(m,m),repmat(Pz,m,1));

%% Solve Bellman equation
 [v,x,pstar] = solvedp(u,P,beta_bah,'policy');  % Bewley-Aiyagari-Huggett
% [v,x,pstar] = solvedp(u,P,beta);                % Uzawa-Epstein
 clear P u C;
%% Steady State Distribution
 d = ergdist(pstar);

%% Summary statistics 
  cTx = YT+R*B-b(x)-AT;  
  pN = (1-a)/a*(cTx./cN).^(1+mu);

 zmean = ergdist(Pz)'*z;
 bmean = b(x)'*d;
 cmean = c'*d;

%% Plot some model properties
 plotdp(v,x,pstar,YT,B,z*yT,b);

%% Simulation
 T      = 500;      
 s0     = findnearest(bmean,B);   
 spath  = simulmarkov(pstar,T,s0);

 yTpath  = YT(spath);
 cTpath  = yTpath + R*B(spath)-b(x(spath))- AT;
 pNpath  = (1-a)/a*(cTpath./cN).^(1+mu);
 CApath = b(x(spath))-B(spath);

 figure(2)
 plot([yTpath cTpath])

 tau =3;
 [mn,sd,corr,acorr] = samplemoms([yTpath cTpath B(spath) CApath 1./pNpath],1,tau);

%% Model steady state statistics  
fprintf('\n\n\n')
fprintf('                     Small Open Economy RBC Model         \n')
fprintf('                  calibrated for the Mexican economy     \n')
fprintf('                      Great Ratios (pct of GDP)           \n')
fprintf('\nRatio              Mexico  Det.SS   Ergodic Simulated  Level')
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

