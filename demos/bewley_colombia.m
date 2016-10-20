%% BAH_UE.M  Bewley-Aiyagari-Huggett vs. Uzawa-Epstein model
%            Follows Durdu, Mendoza and Terrones (2007) NBER-WP 13123
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nDurdu, Mendoza and Terrones (2007) NBER-WP 13123')
 fprintf('\nSingle good SOE model calibrated for Colombia \n')

%% Parameters                       % Mexico calibration values:
 sigma	= 4;                        % risk aversion coefficient 2
 yss    = 0.95;                        % GDP trend mean            1
 R      = 1.035;                    % gross asset return rate   1.059
 phi    = -.4;                      % ad-hoc borrowing limit    -0.51
 
%% Data to match
 bss  = -0.3;                       % NFA to GDP ratio          -0.44
 css  = 0.786;                      % consumption to GDP        0.692
 abs  = yss+bss*(R-1)-css;          % absorption
 rhoy = 0.76;                       % GDP cycle autocorrelation 0.597
 sdvy = 0.026;                      % GDP cycle std deviation   0.03301
  
%% Discount factor
 beta_bah = 0.96;                    % Bewley-Aiyagari-Huggett   0.94
 disc_bah = (1/beta_bah)-1;
 disc_ue  = log(R)/log(1+css);       % Uzawa-Epstein disc rate   
 beta_ue  = (1+css)^-disc_ue;        % Uzawa-Epstein disc factor

%% Stationarity condition
 if R*beta_bah>=1; display('Set beta*R<1 for convergence'); break; end;

%% Markov chain for y
 sdve  = sqrt(sdvy^2*(1-rhoy^2));
 p     = (1+rhoy)/2;
 q     = p;
 ny    = 5;

 [y,Py] = markovchain(ny,p,q,sdvy,yss);

%% State-space S = YxB
 b     = linspace(phi,0.15,1000)';   
 [Y,B] = gridmake(y,b);

 n     = length(y)*length(b); 
 m     = length(b);

%% Utility function and feasible consumption C>0
 c = zeros(n,m);

 for i=1:m    
   c(:,i)=Y+R*B-b(i)-abs;  
 end

 c(c<=0) = NaN;
 u = (c.^(1-sigma)-1)./(1-sigma);
 
 beta = exp(-disc_ue*log(1+c));
 beta(beta>=1) = NaN;

%% Transition probability matrix
 P = kron(speye(m,m),repmat(Py,m,1));

%% Solve Bellman equation
 [v,x,pstar] = solvedp(u,P,beta_bah,'policy');  clear P u c;
% [v,x,pstar] = solvedp(u,P,beta);  clear P u C;

%% Steady State Distribution
 d = ergdist(pstar);
 
%% Summary statistics 
 c = Y+R*B-b(x)-abs;
 
 ymean = ergdist(Py)'*y;
 bmean = b(x)'*d;
 cmean = c'*d;

%% Plot some model properties
 plotdp(v,x,pstar,Y,B,y,b);

%% Simulation
 T    = 5000;      
 s0   = findnearest(bmean,B);   
 s_t  = simulmarkov(pstar,T,s0);

 y_t  = Y(s_t);
 c_t  = y_t + R*B(s_t)-b(x(s_t)) - abs;
 CA_t = b(x(s_t))-B(s_t);

 tau = 3;
 [mn,sd,corr,acorr] = samplemoms([y_t c_t B(s_t) CA_t],1,tau);

%% Model steady state statistics  
fprintf('\n\n\n')
fprintf('                     Small Open Economy RBC Model         \n')
fprintf('                 calibrated for the Colombian economy     \n')
fprintf('                      Great Ratios (pct of GDP)           \n')
fprintf('\nRatio             Colombia  Det.SS   Ergodic Simulated  Level')
fprintf('\nOutput           %8.2f %8.2f %8.2f %8.2f %8.2f'   ,1,yss/yss,ymean/ymean,mn(1)/mn(1),mn(1))
fprintf('\nConsumption      %8.2f %8.2f %8.2f %8.2f %8.2f'   ,0.786,css/yss,cmean/ymean,mn(2)/mn(1),mn(2))
fprintf('\nNFA to output    %8.2f %8.2f %8.2f %8.2f %8.2f'   ,bss  ,bss/yss,bmean/ymean,mn(3)/mn(1),mn(3))
fprintf('\nTime b.c.        %8.2f %8.2f %8.2f %8.2f %8.2f\n' ,0.162,0,sum(d(1:ny)),NaN,NaN)

fprintf('\n')
fprintf('                                  Moments                    \n')
fprintf('                            Data                        Model  ')
fprintf('\nVariable          Std.Dev  Corr(i,y) Autocorr  Std.Dev Corr(i,y) Autocorr')
fprintf('\nOutput           %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,2.47,1.00,0.76,sd(1),corr(1,tau+1),acorr(tau,1))
fprintf('\nConsumption      %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,3.31,0.83,0.80,sd(2),corr(2,tau+1),acorr(tau,2))
fprintf('\nNFA              %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'   ,NaN ,NaN ,NaN ,sd(3),corr(3,tau+1),acorr(tau,3))
fprintf('\nCurrent account  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' ,NaN ,NaN ,NaN ,sd(4),corr(4,tau+1),acorr(tau,4))

