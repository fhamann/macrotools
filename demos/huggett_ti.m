% TIME ITERATION HUGGETT TOOLBOX
  disp('PRUEBA 1 TIME ITERATION')
  clear all  

% ENTER MODEL PARAMETERS
beta    = 0.993;    % Discount factor
gamma   = 1.5;        % Risk aversion
blim    = -8;       % Ad hoc debt limit
q       = 0.997;    % exogenous price of debt 

% SET THE EXOGENOUS STATES AND THE TRANSITION PROB. MATRIX
ne      = 2;        % Number of exogenous prod. states
[e,Pe]  = markovchain(ne,0.5,0.8,0.1,0.5);

% PACK MODEL STRUCTURE
  clear model
%  model.func = 'hug_func';               % model function file
  model.func = 'hug_func_inv';           % model function file w. inverse
  model.discount = beta;                 % discount factor
  model.e = e;                           % shocks
  model.w = Pe;                          % probabilities
  model.params = {gamma q};              % other parameters
  
% DEFINE APPROXIMATION SPACE
  n      = 50;                             % degree of approximation
  bmin   = blim;                          % minimum state
  bmax   = -blim*3/4;                     % maximum state
  fspace = fundefn('lin',n,bmin,bmax);   % function space
  bnodes = funnode(fspace);              % state collocaton nodes

% SOLVE MODEL

  [c,bp,v,resid] = solveti(model,fspace,bnodes);


figure(1)
plot(bnodes,bp,bnodes,bnodes,'--')
legend('policy in state 1','policy in state 2','$45^o$ line','Location','southeast')
title(['Huggett Partial Equilibrium Model: $\gamma=$', num2str(gamma),...
       ', $\beta=$',num2str(beta),', Q =', num2str(q)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Net Assets $a_{t}$','interpreter','latex','FontSize',12)
ylabel('Net Assets choice $a_{t+1}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
 
    figure(2)
   plot(bnodes,v);
legend('Value in state 1','Value in state 2','Location','southeast')
title(['Huggett Partial Equilibrium Model: $\gamma=$', num2str(gamma),...
       ', $\beta=$',num2str(beta),', Q =', num2str(q)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Net Assets $a_{t}$','interpreter','latex','FontSize',12)
ylabel('Value $v_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
 
%% SIMULATE MODEL

[hs,hx,hs2x] = simulti(e,bnodes,bp,c,Pe);
  