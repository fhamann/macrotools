function [amean,phi,d]=Ea(rate,beta,mu,delta,alpha,b,tau,zmean,rho,sigma)

% Ea.m Aiyagari's asset demand function under incomplete markets

% Set number of actions and states 
 nz = 7;                                 % # gridpoints exog. states
 na = 100;                               % # gridpoints endo. states
 n  = nz*na;                             % # total states gridpoints
 m  = na;                                % # total action gridpoints

%% Construct State-Space
% 1. Approximate shocks with n1 states Markov chain
 [z, prob] = rouwenhorst(nz,zmean,rho,sigma);
 z = exp(z);

% 2. Asset Grid (and setting the borrowing limit)
 wage = (1-alpha)*(zmean*(alpha/(rate+delta))^alpha)^(1/(1-alpha));

 if rate<=0
    phi = b;                             % ad-hoc borrowing limit
 else
    phi  = min(b, wage*z(1)/rate);       % natural borrowing limit
 end

 a = linspace(-phi,40,na);               % linearly spaced vector

 [Z,A] = gridmake(z',a');

%% Construct the reward function, u
 C  = zeros(n,m);

  for i=1:m
    C(:,i) = (wage*Z+(1+rate)*A-a(i)-tau);
  end

 C(C<=0)=nan;
 u = (C.^(1-mu)-1)/(1-mu); 

%% Construct the transition matrix
 P  = kron(speye(na,na),repmat(prob,na,1));                   

%% Solve the Model
 [v,x,pstar] = solvedp(u,P,beta,'policy');

%% Ergodic Distribution
  d = ergdist(pstar);
  amean =a(x)*d;