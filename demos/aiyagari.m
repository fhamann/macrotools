% AIYAGARI Aiyagari (1993) model
%
% Written by F. Hamann. Feel free to copy, change and distribute
 fprintf('\nAiyagari 1993 model \n')

%% Model Parameters
 mu        = 3;                         % risk aversion              
 discrate  = 0.05;                      % discount rate
 beta      = 1/(1+discrate);            % subjective discount factor 
 delta     = 0.08;                      % depreciation
 Z         = 1;                         % production technology
 alpha     = 0.36;                      % capital's share of income
 b         = 3;                         % ad-hoc debt limit
 zmean     = 1;                         % mean of shocks
 rho       = 0.2;                       % autocorrelation coefficient
 sigma     = 0.4*sqrt(1-rho*rho);       % standard deviation of shocks
 tau       = 0.1;                       % tax rate

%% Markov chain 
 [z, prob] = rouwenhorst(7,zmean,rho,sigma);        
 z = exp(z);

%% Plot r vs. Ea (equilibrium assets)
 r      = [0:0.0025:discrate-0.0025];   
 amean  = zeros(1,length(r));
 phi    = zeros(1,length(r));

 for i=1:length(r)
  [amean(i),phi(i)] = Ea(r(i),beta,mu,delta,alpha,b,tau,zmean,rho,sigma);
 end

 figure(1)
 plot(amean,r); hold on;
 plot(tau./r,r); hold on;
 plot(-phi,r)
 title('Aiyagari model graphical interpretation')
 
%% Equilibrium r y Ea

 % Bisection method:
 % 1.) Define initial r, r0, and compute Ea(r0)
 % 2.) Find r1 such that bhat(r1)=Ea(r0)
 % 3.) Define r2=(r0+r1)/2 and compute Ea(r2) (bisection)
 % 4.) If Ea(r2)>bhat(r2), then replace r0 by r2;
 %     if Ea(r2)<bhat(r2), then replace r1 by r2,
 %     and use bisection.

 r0=0.04 % r0 < discount rate, but close
 [am0,phi0]=Ea(r0,beta,mu,delta,alpha,b,tau,zmean,rho,sigma);
 r1=tau/am0;

 maxit = 100; tol = 10e-4;
 
 for i=1:maxit
    r2=(r0+r1)/2;
    [am2,phi2] = Ea(r2,beta,mu,delta,alpha,b,tau,zmean,rho,sigma);
    bhat2=tau/r2;i,am2,bhat2
    if abs(am2-bhat2)<tol,break,end
    if am2>bhat2
        r0=r2;
    else am2<bhat2
        r1=r2;
    end 
 end    

 disp('Equilibrium interest rate')  , r2
 disp('Equilibrium stock of assets'), am2

%% Stationary distribution of assets at equilibrium rate
 [Eae,phi,d] = Ea(r2,beta,mu,delta,alpha,b,tau,zmean,rho,sigma);

 hold off
 figure(2); bar(d)