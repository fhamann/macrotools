%% OCCHOICE.M  Occupational choice model
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nOccupational choice model \n')

%% Parameters
 sigma	= 2;            % risk aversion
 beta	= 0.99;         % discount factor
 R      = 1.004;        % gross asset return rate
 tau    = 0.2;          % income tax
 kappa  = 0.1;          % switching cost

 if R*beta>=1; display('Set beta*R<1 for convergence'); break; end;

%% Markov chain for y: [y,prob] = markovchain(ny,p,q,e,m,a)
 [y,Py] = markovchain(3,0.7,0.7,0.6);      % Markov chain
 y      = exp(y);                          % set earnings positive
 yf     = y'*ergdist(Py);                  % formal earnings

%% State-space S = YxB
 b       = linspace(-15,10,300)';          % B = {b1<b2<b3<...<bn}
 f       = [0;1];                          % F = {0,1}; 0:formal,1:informal

 [Y,B,F] = gridmake(y,b,f);                % state-space grid
 [bp,fp] = gridmake(b,f);                  % action-space grid

 n       = length(y)*length(b)*length(f);  % number state gridpoints
 m       = length(b)*length(f);            % number action gridpoints  

%% Utility function and feasible consumption C>=0
 C  = zeros(n,m);
  for i=1:m
      C(:,i) = ((1-tau)*yf + R*B - bp(i)).*(1-F)*(1-fp(i)) + ...
               ((1-tau)*yf + R*B - bp(i)).*(1-F)*fp(i)     + ...
               (Y+R*B-bp(i)).*F*fp(i)                      + ...
               (Y+R*B-bp(i)-kappa).*F*(1-fp(i));
  end

 C(C<=0) = NaN;
 u  = (C.^(1-sigma)-1)./(1-sigma);

%% Transition probability matrix (see Sargent and Ljundqvist)
 P  = kron(speye(m,m),repmat(Py,m,1));      

%% Bellman equation
 [v,x,pstar] = solvedp(u,P,beta,'policy'); clear P u C;

%% Steady State Distribution
 d = ergdist(pstar); plot(d)

%% Print results  
  fprintf('\nSteady-state Model Statistics \n ')
  fprintf('\nOccupational choices ')
  fprintf('\n Formal             %8.2f'  ,sum(d(1:n/2)))  
  fprintf('\n Informal           %8.2f'  ,sum(d(1+n/2:end))) 
  fprintf('\nNet Assets          %8.2f'  ,bp(x)'*d) 
  fprintf('\n Formal             %8.2f'  ,bp(x(1:n/2))'*d(1:n/2))  
  fprintf('\n Informal           %8.2f\n',bp(x(1+n/2:end))'*d(1+n/2:end)) 