% EATON_GERSOVITZ Hamann (2002) model by discretization
%
% This code solves the E&G model for trend shocks and exogenous labor
% supply. The notation and some assumptions here differs slightly from 
% those on the paper.

clear,clc

% Model Parameters
mu     = 2;                              % risk aversion              
beta   = 0.8;                            % subjective discount factor 
delta  = 0.02;                           % output loss in autarky
alpha  = 0.32;                           % capital's share of income
lambda = 0.1;                            % prob of redemption 0.1
rate   = 0.01;                           % exogenous interest rate

n1 = 15;                                  % # gridpts exog. for growth rate 
n2 = 200;                                % # gridpts endo. for assets (400)
n  = n1*n2;                              % # total states gridpoints
m  = n2;                                 % # total action gridpoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the State-Space                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Approximate Ln(g_t)= (1-rho)*(Ln(meg)-cg) + rho*Ln(g_t-1) + e_t
%    with n1 states Markov chain. Set width of state shocks larger
%    than because default ocurrs in extreme states

rho    = 0.17;                           % autocorrelation coefficient
sigma  = 0.03;                           % std dev of growth rate shocks
meg    = 1.006;                          % long run mean for trend income
cg     = .5*((sigma^2)/(1-rho^2));       % AR(1) parameter
cg     = (1-rho)*(log(meg)-cg);          % constant term of the AR(1)
width  = 4.1458;                         % width of state shocks
                                         
[prob,g,pi] = appmarkov2(cg,lambda,sigma,width,n1);
g = exp(g);

% 2. Build the state-space grid 

amax  = 0;                                % max value of asset grid  
amin  = -0.22;                            % min value of asset grid  
a     = linspace(amin,amax,n2);           % linearly spaced vector
azero = find(a==0);                       % state where assets is zero

[A,G] = gridmake(a',g');                 % the order of states follows A&G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer-loop: find "q" that satisfies the capital market equilibrium      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Initial guess for q

q0 = (1/(1+rate))*ones(n,m);

maxit = 1000;
for it = 1:maxit
  q0old = q0;    

  % 1. Construct the reward function without default and solve

  y  = G./meg;
  y  = repmat(y,1,m);  

  c  = zeros(n,m);
  f  = zeros(n,m);

  for i=1:m
    c(:,i) = y(:,i)+A-(G*a(i)).*q0(:,i); 
  end

  f = (c.^(1-mu))/(1-mu); f(c<=0) = NaN;

  % 2. Construct the reward function with default 

  cdefault  = zeros(n,m);
  fdefault  = zeros(n,m);

  for i=1:m
    cdefault(:,i) = (1-delta)*y(:,i);
  end

  fdefault = (cdefault.^(1-mu))/(1-mu); fdefault(cdefault<=0) = NaN;

  % 3. Comupte new discount factor to account for the detrend form

  betah = (G.^(1-mu))*beta;
  betah = repmat(betah,1,m);

  % 4. Solve the system of Bellman eq

  [vd,xd,vr,xr,D] = solvedpos(f,fdefault,prob,betah,lambda,azero);

  % 5. Compute the expected value of default

  Edef = prob*(reshape(D,m,n1)');
  ind  = reshape(D,m,n1)';
  ind  = sum(ind)==n1;                 % find a' where default is Prob=1 

  for j=1:m;
      if ind(j)==1;
        Edef(:,j)=1;
      end;
  end;

  Edef = kron(Edef,ones(m,1)); % add a(t) to rows to make same size as q

  % 6. Update the price "q"

  q0=(1/(1+rate))*(1-Edef);
  q0=max(q0,0);

  if abs(q0old-q0)<10e-6,break,end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of iteration to find the "q"                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10. Optimal Policy Function

% x = (1-D).*xr+D.*azero;

x     = (1-D).*xr+D.*xd;          % policy function    
Pstar = otpm(x,prob,n2,1,n1);     % optimal transition prob matrix
pie   = ergdist(Pstar);           % ergodic distribution

Psr = otpm(xr,prob,n2,1,n1);
Psd = otpm(xd,prob,n2,1,n1);

Eq    = (pie'*q0)';
Ebp   = pie'*A(x);
Dr    = pie'*D(x);

figure(1);
subplot(1,3,1); plot(pie); 
subplot(1,3,2); spy(Pstar); 
subplot(1,3,3); plot(a,Eq)

%% Tablas

fprintf('\n\n\n      Sovereign Default Model \n') 
fprintf('\n                Parameters        \n') 
fprintf('\nDiscount factor                    %5.2f'   ,beta)
fprintf('\nRisk-free rate                     %5.2f'   ,rate*100)
fprintf('\nRisk aversion                      %5.2f'   ,mu)
fprintf('\nProbability of redemption          %5.2f'   ,lambda)
fprintf('\nDefault penalty                    %5.2f'   ,delta*100)
fprintf('\nBorrowing limit (pct)              %5.1f'   ,amin*100)
fprintf('\nEndowment std dev (pct)            %5.2f'   ,sigma*100)
fprintf('\nEndowment autocorrelation          %5.2f'   ,rho)
fprintf('\n')
fprintf('\n           Ergodic Moments          \n') 
fprintf('\nNumber of iterations over q        %5.0f'   ,it)
fprintf('\nExpected Net Foreign Assets        %5.2f'   ,Ebp)
fprintf('\nErgodic default rate (pct)         %5.2f'   ,Dr*100)
fprintf('\n')