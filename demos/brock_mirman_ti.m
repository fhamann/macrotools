% TIME ITERATION BM TOOLBOX
  disp('PRUEBA 2 TIME ITERATION')
  clear all  

% ENTER MODEL PARAMETERS
beta    = 0.98;    % Discount factor
alpha   = 0.3;

% COMPUTE GAUSSIAN NODES AND WEIGHTS
 ne      = 3;        % Number of exogenous prod. states
 zmean   = 1;       
 zstdv   = .2;      
 [e,Pe]  = markovchain(ne,0.5,0.5,zstdv,zmean);

% PACK MODEL STRUCTURE
  clear model
model.func = 'brock_mirman_func';      % model function file
%   model.func = 'brock_mirman_func_inv';  % model function file w. inverse
  model.discount = beta;                 % discount factor
  model.e = e;                           % shocks
  model.w = Pe;                          % probabilities
  model.params = {alpha};                % other parameters
  
% DEFINE APPROXIMATION SPACE
  n      = 50;                              % Number of grid points for the state variable
  kstar  = (alpha*beta)^(1/(1-alpha));      % SS for the State variable
  kmin   = kstar*0.2;
  kmax   = kstar*1.8;
  fspace = fundefn('lin',n,kmin,kmax);
  k      = funnode(fspace);  

% EXACT SOLUTION
   [K,E] = gridmake(k,e);
   kex = reshape((alpha*beta*E.*K.^alpha),n,ne);

  [x,kp,v] = solveti(model,fspace,k);

figure(1)
plot(k,kp,k,kex,'x')
legend('policy in state 1','policy in state 2','policy in state 2','exact solution','Location','southeast')
title(['Huggett Partial Equilibrium Model: $\alpha=$', num2str(alpha),...
       ', $\beta=$',num2str(beta)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Capital $k_{t}$','interpreter','latex','FontSize',12)
ylabel('Capital choice $k_{t+1}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')

%% SIMULATE MODEL

[hs,hx,hs2x] = simulti(e,k,kp,x,Pe);