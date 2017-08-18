% AGUIAR-GOPINATH Aguiar & Gopinath (2006) Model II by state space discretization
% The present code solves the A&G model for a stochastic trend and exogenous labor
% supply. The notation here is not the same used in the A&G code.

clear
clc
tic
% Model Parameters
mu     = 2;                                    % risk aversion              
beta   = 0.8;                                  % subjective discount factor 
delta  = 0.02;                                 % additional lost of output in autarky
alpha  = 0.32;                                 % capital's share of income
lambda = 0.1;                                  % prob of redemption 0.1
A      = 1;                                    % Output scale
rate   = 0.01;                                 % exogenous interest rate

% Parameters for the exogenous state variables grids
rho_z    = 0.9;                                % autocorrelation coefficient
sigma_z  = 0.0;                                % standard deviation of shocks for growth rate
mn_z     = (-1/2)*(sigma_z)^2;                 % long run mean for trend income
rho_g    = 0.17;                               % autocorrelation coefficient                                 0.17
sigma_g  = 0.03;                               % standard deviation of shocks for growth rate                0.03
mn_g     = 1.006;                              % long run mean for trend income                             1.006
cons_g   = (1-rho_g)*(log(mn_g)-(.5*((sigma_g^2)/(1-rho_g^2))));                % constant term of the AR(1)
width    = 4.1458;                             % width of state shocks, needs to be larger becasue default ocurrs in extreme states
n1       = 25;                                 % # gridpoints exo. states for growth

% Parameters for the endogenous state variable grid
bmax   = 0;                                    % max value of asset grid  
bmin   = -0.22;                                % min value of asset grid  
n2     = 400;                                  % # gridpoints endo. states for assets (400)

% Computational parameters
n      = n1*n2;                                % # total states gridpoints
m      = n2;                                   % # total action gridpoints
maxit  = 10000;                                % maximum iteration number

%% Construct State-Space
% Exogenous state variable 
% 1. Approximate the trend shock as:
%   Ln(g_t)= (1-rho)*(Ln(meg)-cg) + rho*Ln(g_t-1) + e_t
% with n1 states Markov chain. The first model is a stable trend,
% sigma_g = 0; rho_g = NA;

if sigma_g >0
[gprob,g] = appmarkov2(cons_g,rho_g,sigma_g,width,n1); g = g';
trend = exp(g);
else 
trend = mn_g;     
end

% 2. Approximate the transitory component
if sigma_z>0
[z,zprob] = tauchenhussey(n1,mn_z,rho_z,sigma_z,sigma_z);
% z = linspace(mn_z-2.5*sigma_z,mn_z+2.5*sigma_z,n1)';
else
z = 0;
end
y = A*exp(z)*trend;

% 3. Asset Grid 

b     = linspace(bmin,bmax,n2)';               % linearly spaced vector
incb  = (max(b)-min(b))/(n2-1);                % width of the asset grid
bzero = find(b==0);                            % state where assets is zero

% 4. Full State space

[B,Z,TR] = gridmake(b,z,trend);                % the order of states follows A&G code
Y = A*exp(Z).*TR;

% 5. Variable Detrending

Y_bar = Y./(mn_g*TR);
B_bar = B./(mn_g*TR);
b_bar = b/mn_g;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration to find the "q" that satisfies the capital market equilibrium %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6. Initial guess for q

q = (1/(1+rate))*ones(n,m);

cpay  = zeros(n,m);         cdef  = zeros(n,m);
upay  = zeros(n,m);         udef  = zeros(n,m);

fprintf('\nAGUIAR-GOPINATH STOCHASTIC TREND MODEL     \n') 
for it = 1:maxit
    
    qold = q;    

% 7. Construct the reward function with default (def) and without default (pay) 

  for ib = 1:m
    cpay(:,ib) = Y_bar + B_bar - q(:,ib).*b_bar(ib); 
    cdef(:,ib) = (1-delta)*Y_bar;
  end

upay = (cpay.^(1-mu))/(1-mu); upay(find(cpay<=0)) = NaN;
udef = (cdef.^(1-mu))/(1-mu); udef(find(cdef<=0)) = NaN;

% 8. Solve the Dynamic Programming problem
   fprintf('Iteration %2.0f', it)
[vdef,xdef,vpay,xpay,V,D] = solvedpos(upay,udef,gprob,beta,lambda,bzero);

% 9. Compute the expected value of default

Edef = gprob*(reshape(D,m,n1)');
ind = reshape(D,m,n1)';
ind = sum(ind)==n1; %ind finds a' in which default happens with prob 1 next period

for j=1:m;
    if ind(j)==1;
        Edef(:,j)=1;
    end;
end;

Edef = kron(Edef,ones(m,1)); %add current a to rows to make same size as q

% 10. Update the price "q"

q = (1/(1+rate))*(1-Edef);   q = max(q,0);

if abs(qold-q)<10e-6,break,end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of iteration to find the "q"                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optimal policy & transition prob matrix, for the equilibrium "q"
x       = (1-D).*xpay+D.*bzero;                   % policy function at b'=0 
Pstar   = otpm1v(x,gprob,n2,n1);
pie     = ergdist(Pstar);                            % ergodic distribution
Eq      = (pie'*q)';                                              
Ebp_bar = pie'*b_bar(x);
Ebp     = pie'*b(x);

% Before simulating the model, define an approximation space for the bond price
qspace  = fundefn('lin',35,bmin,bmax);
[cq,Bq] = funfitxy(qspace,b,Eq);

simul_aguigop