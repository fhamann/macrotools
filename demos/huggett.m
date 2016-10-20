% HUGGETT.M  PARTIAL EQUILIBRIUM version of Huggett (1993)
%            See Huggett (1993) The importance of the risk-free rate
%            in heterogeneous-agent incomplete insurance economies
%            Journal of Economic Dynamics and Control 17 953-969
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nHeterogeneuos agents and incomplete markets model \n')

%% Parameters
 sigma = 1.5;
 beta  = 0.993;
 q     = 0.997;

 if beta>=q; display('Set beta<q for convergence'); break; end;

%% Markov chain for e: [e,Pe] = markovchain(ny,p,q,eps,m,a)
% [e,Pe] = markovchain(2,0.5,0.8,0.5,0.5,0.8);      % Markov chain
 [e,Pe] = markovchain(2,0.5,0.8,0.5,1);      % Markov chain
 
%% State-space S = ExA  
 a = linspace(-8,3,250)';

 n = length(e)*length(a); 
 m = length(a);

 [E,A] = gridmake(e,a);
 
%% Utility function and feasible consumption C>=0
 C = zeros(n,m);
 
 for i=1:m
    C(:,i) = E+A-q*a(i);
 end
 C(C<=0) = nan;

 u = (C.^(1-sigma)-1)./(1-sigma);
 
%% Transition probability matrix (see Sargent and Ljundqvist)
 P = kron(speye(m,m),repmat(Pe,m,1));

%% Bellman equation
 [v,x,pstar] = solvedp(u,P,beta,'policy'); clear P u C;

%% Steady State Distribution
 d  = markov(pstar);
 c = E+A-q.*a(x);

%% Summary statistics
 amean = a(x)'*d;
 cmean = c'*d;
 emean = E'*d;
 sd.c  = sqrt(((E+A-q*a(x)-cmean).^2)'*d);
 sd.e  = sqrt(((E-emean).^2)'*d);

%% Plot figures
 plotdp(v,x,pstar,E,A,e,a)

%% Print results  
 fprintf('\nSteady-state Model Statistics \n ')
 fprintf('\nPopulation means ')
 fprintf('\n Earnings           %8.2f'  ,emean)  
 fprintf('\n Consumption        %8.2f'  ,cmean) 
 fprintf('\n Net Assets         %8.2f'  ,amean) 
 fprintf('\n Assets to income   %8.2f'  ,amean/emean) 
 fprintf('\nPopulation volatility') 
 fprintf('\n Consumption        %8.2f'  ,sqrt(((E+A-q*a(x)-cmean).^2)'*d)) 
 fprintf('\n Earnings           %8.2f\n',sqrt(((E-emean).^2)'*d))