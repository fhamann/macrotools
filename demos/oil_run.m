%% OPTIMAL OIL EXTRACTION WITH FIXED DISCOVERIES

clear all;

%% MODEL PARAMETERS
   
beta  = 1/1.04;           % discount factor                             1/R
phi   = 4.25;             % extraction cost parameter                  4.25
gama  = 0.99;             % extraction cost parameter                  0.99
d     = 0.056;            % fixed discovery rate                      0.056

%% TECHNICAL PARAMETERS

T        = 100000;        % # of periods of simulation

%% LOAD DATA AND SET OIL PRICE SPACE

oil_price_ts;                                     % load WTI data -> "poil"

lp        = log(poil);
lpt       = hpfilter(lp,1600);                
lpc       = lp-lpt;
[rho,sig] = ols(lpc(1:end-1),lpc(2:end));
[logp,Q]  = tauchenhussey(2,mean(lpc),rho,sig,sig);
dp        = markov(Q);
pnodes    = (exp(logp')./(exp(logp')*dp))';

%% PACK MODEL STRUCTURE

%model.func     = 'oilcomp_func';
model.func     = 'oilcomp_func_inv';
model.discount = beta;
model.e        = pnodes;
model.w        = Q;
model.params   = {phi gama d};

%% DEFINE APPROXIMATION SPACE

n      = 15;                           % degree of approximation
smin   = 0;                            % minimum state
smax   = 1;                            % maximum state
fspace = fundefn('lin',n,smin,smax);   % function space
snodes = funnode(fspace);              % state collocaton nodes

%% SOLVE MODEL

[x,sp,v,resid] = solveti(model,fspace,snodes);

%% PLOT FIGURES

figure(1)
plot(snodes,sp,snodes,snodes,'--')
legend('policy in state 1','policy in state 2','$45^o$ line','Location','southeast')
title(['Optimal Oil Extraction Model: $\gamma=$', num2str(gama),...
       ', $\beta=$',num2str(beta),', $\phi=$', num2str(phi)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Oil reserves $s_{t}$','interpreter','latex','FontSize',12)
ylabel('Next periods oil reserves $s_{t+1}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
 
figure(2)
plot(snodes,v);
legend('Value in state 1','Value in state 2','Location','southeast')
title(['Optimal Oil Extraction Model: $\gamma=$', num2str(gama),...
       ', $\beta=$',num2str(beta),', $\phi=$', num2str(phi)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Oil reserves $s_{t}$','interpreter','latex','FontSize',12)
ylabel('Value $v_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')

%% SIMULATE MODEL

[hs,hx,hs2x] = simulti(pnodes,snodes,sp,x,Q);
