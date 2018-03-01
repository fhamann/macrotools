% EATON_GERSOVITZ Hamann (2002) implementation by discretization
%
% This code solves the E&G model for trend shocks and exogenous labor
% supply. The notation and some assumptions here differs slightly from 
% those on the paper.

clear

%% Model parameters 

% Parameters set to match model moments
beta    = 0.825;                  % Subjective base discount factor
lambda  = 1/(2*4);                % Re-entry prob (lambda = 1/(autyears*4)) Colombia 3
d0      = 0.025;                  % Level parameter for default cost
d1      = 1.5;                    % Curvature parameter for default cost

% "Free" or conventional parameters
mu      = 3;                      % risk aversion              
rf      = 0.01;                   % exogenous interest rate

% Values computed from Colombian data
mn_g    = 1.008613935676046;      % long run mean for trend income 

% Taking into account that z = e+m
rho_z   = 0.97;                   % autocorrelation coefficient for the persistent part
sigma_z = 0.000091469;            % var of productivity shocks persistent part
% sigma_m = 0.000025;               % var of productivity shocks purely transitory part

embi    = 0.03482;                % 0.03482 
c2GDP   = 0.8451;                 % 0.8451
a2GDP   = 0.2907;                 % 0.2907 Total debt // 0.1732 Public Debt

% Data implied parameters
r       = rf+embi;                % Colombian interest rate
X       = 1+r*a2GDP-c2GDP;        % Investment (constant)

% State space parameters
amax    = 0;                      % max value of asset grid  
amin    = -0.6;                   % min value of asset grid  
ng      = 1;                      % # gridpts for growth rate 
nz      = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the State-Space                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Build the exogenous state-space 

    g      = mn_g;
    prob_g = 1;

[prob_z,z] = appmarkov2(0,rho_z,sqrt(sigma_z),3,nz);  z = exp(z);  

prob_gz = kron(prob_g,prob_z);
prob_gz = prob_gz./(((sum(prob_gz,2)))*ones(1,nz*ng));

% 2. Build the endogenous state-space grid 

% a     = amin:spa:amax;          % Asset grid  
a = linspace(amin,amax,450);
azero = find(a==0);             % state where assets is zero

% 3. State-space
[A,Z,G] = gridmake(a',z',g');     

nex = nz*ng;                    % Exogenous states number of gridpts
na  = length(a);                % Endogenous states number of gridpts
n   = na*nex;                   % State number of nodes
m   = na;                       % Actions number of nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer-loop: find "q" that satisfies the capital market equilibrium      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y   = (Z.*G)./mn_g;
y   = repmat(y,1,m);  
cd  = (1-(d0*repmat(Z,1,m).^d1)).*y-X;       % Consumption under default
ud  = (cd.^(1-mu))/(1-mu); ud(cd<=0) = NaN;  % Utility under default 

cnd = zeros(n,m);   
und = zeros(n,m);

% 0. Initial guess for q

q = (1/(1+rf))*ones(n,m);

maxit = 1000;
for it = 1:maxit
  qold = q;    

  % 1. Construct the reward function with and without default

  for i=1:m
    cnd(:,i) = y(:,i)+A-(G*a(i)).*q(:,i)-X;       % Consumption under no default
  end

  und = (cnd.^(1-mu))/(1-mu); und(cnd<=0) = NaN; % Utility under no default

  % 3. Comupte new discount factor to account for the detrend form

  beta_h = (G.^(1-mu))*beta;
  beta_h = repmat(beta_h,1,m);

  % 4. Solve the system of Bellman eq

  [vd,xd,vr,xr,D,V] = solvesovereign(und,ud,prob_gz,beta_h,lambda,azero);
  
  % 5. Compute the expected value of default

  Edef = prob_gz*(reshape(D,m,nex)');
  ind  = reshape(D,m,nex)';
  ind  = sum(ind)==nex;                  

  for j=1:m;
      if ind(j)==1;
        Edef(:,j)=1;
      end;
  end;

  Edef = kron(Edef,ones(m,1)); % add a(t) to rows to make same size as q

  % 6. Update the price "q"

  q=(1/(1+rf))*(1-Edef);
  q=max(q,0);

  fprintf(' Change in q = %5.6f  after iteration %5.0f\n',norm(q-qold),it) 

  if abs(qold-q)<10e-6,break,end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of iteration to find the "q"                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10. Optimal Policy Function

x = (1-D).*xr+D.*azero;

Pstar = otpm(x,prob_gz,na,1,nex);     % optimal transition prob matrix
pie   = ergdist(Pstar);           % ergodic distribution

Eq    = (pie'*q)';
Eap   = pie'*A(x);
Dr    = pie'*D;

for i=1:n;
     cp(i,1) = cnd(i,x(i));
     qp(i,1) = q(i,x(i));
end;

% %% Simulation 
% 
    N = 1;      % number of simulations
    T = 250000; % Sample simulation
    J = 0;      % Drop first J simulations
    L = 1;      % Dynamic lags for cross-correlations
%    
    s0 = getindex([Eap mn_g],[A G]);
% 
    nv = 11;   % Number of variables to compute sample moments
%      
    M  = zeros(nv,1,N);      Md  = zeros(nv,1,N);     Mp  = zeros(nv,1,N);
    SD = zeros(nv,1,N);      SDd = zeros(nv,1,N);     SDp = zeros(nv,1,N);
    XC = zeros(nv,1+2*L,N);  XCd = zeros(nv,1+2*L,N); XCp = zeros(nv,1+2*L,N);
    AC = zeros(1+2*L,nv,N);  ACd = zeros(1+2*L,nv,N); ACp = zeros(1+2*L,nv,N);

    for i = 1:N
   % Initial values 
   At      = zeros(T,1);
   at      = zeros(T,1);
   dt      = zeros(T,1);
   qt      = zeros(T,1);
   ct      = zeros(T,1);
   yt      = zeros(T,1);
   gt      = zeros(T,1);
   zt      = zeros(T,1);
   CAt     = zeros(T,1);
   TBt     = zeros(T,1);
   drst    = zeros(T,1);
   history = zeros(T,1);
   dt      = zeros(T,1);
   
   ia = max(1,round(rand*na)); At(1) = a(ia);
   ig = simulmarkov(Pstar,T,2)';
                
     redemp = ones(floor(lambda*T),1);
     redemp = [redemp; zeros(T-length(redemp),1)];
     temp   = randperm(T)';
     redemp = redemp(temp);

 for t = 1: T   
     if history(t)==0 | (history(t)==1 & redemp(t)==1); 
        dt(t) = D(ig(t));
        
     if dt(t) == 0; 
        
        at(t)        = a(x(ig(t)));
        qt(t)        = q(ig(t),x(ig(t)));
        ct(t)        = cnd(ig(t),x(ig(t)));
        yt(t)        = (Z(ig(t))*G(ig(t)))./mn_g;
        gt(t)        = G(ig(t));
        zt(t)        = Z(ig(t));
        CAt(t)       = A(x(ig(t)))-A(ig(t));
        TBt(t)       = yt(t) - ct(t)- X;
        drst(t)      = ((1./qt(t)-1)-(rf));
                
     elseif dt(t) == 1;  
        at(t)        = 0;
        qt(t)        = 0;
        ct(t)        = cd(ig(t),x(ig(t)));
        yt(t)        = (Z(ig(t))*G(ig(t)))./mn_g;
        gt(t)        = G(ig(t));
        zt(t)        = Z(ig(t));      
        CAt(t)       = A(x(ig(t)))-A(ig(t));
        TBt(t)       = yt(t) - ct(t) - X;
        drst(t)      = 0;
        history(t+1) = 1;
     end
     
 elseif history(t)==1 & redemp(t)==0;
        at(t)        = 0;
        qt(t)        = 0;
        ct(t)        = cd(ig(t),x(ig(t)));
        yt(t)        = (Z(ig(t))*G(ig(t)))./mn_g;
        gt(t)        = G(ig(t));
        zt(t)        = Z(ig(t));
        CAt(t)       = A(x(ig(t)))-A(ig(t));
        TBt(t)       = yt(t) - ct(t) - X;
        drst(t)      = 0;
        history(t+1) = 1;
  end
        At(t+1) = at(t);
 end

 if length(history) > T
     history = history(1:end-1);
 end
 defaulteps = (find((history==1)&(redemp==0)|(history==0)&(dt==1)));
 At = At(1:end-1); 
%  ctt = hpfilter(log(ct),1600); ctc = ct - ctt;
%  ytt = hpfilter(log(yt),1600); ytc = yt - ytt;
%  qtt = hpfilter(qt,1600);      qtc = qt - qtt;
 
 
        [M(:,:,i),SD(:,:,i),XC(:,:,i),AC(:,:,i)]     = samplemoms([At at qt ct dt yt history redemp CAt./yt TBt./yt drst ],6,L);
        [Mp(:,:,i),SDp(:,:,i),XCp(:,:,i),ACp(:,:,i)] = samplemoms([At(find((history==0)&(dt==0))) at(find((history==0)&(dt==0))) qt((find((history==0)&(dt==0)))) ct((find((history==0)&(dt==0)))) dt((find((history==0)&(dt==0)))) yt((find((history==0)&(dt==0)))) history((find((history==0)&(dt==0)))) redemp((find((history==0)&(dt==0)))) CAt((find((history==0)&(dt==0))))./yt((find((history==0)&(dt==0)))) TBt((find((history==0)&(dt==0))))./yt((find((history==0)&(dt==0)))) drst((find((history==0)&(dt==0)))) ],6,L);
        [Md(:,:,i),SDd(:,:,i),XCd(:,:,i),ACd(:,:,i)] = samplemoms([At(defaulteps) at(defaulteps) qt(defaulteps) ct(defaulteps) dt(defaulteps) yt(defaulteps) history(defaulteps) redemp(defaulteps) CAt(defaulteps)./yt(defaulteps) TBt(defaulteps)./yt(defaulteps) drst(defaulteps) ],6,L);
 
    end
%    % Take averages of the N-simulations over the third dimensions
if N>1
    mn    = mean(M,3) ;    mnd    = mean(Md,3) ;    mnp    = mean(Mp,3) ;
    sd    = mean(SD,3);    sdd    = mean(SDd,3);    sdp    = mean(SDp,3);
    xcorr = mean(XC,3);    xcorrd = mean(XCd,3);    xcorrp = mean(XCp,3);
    acorr = mean(AC,3);    acorrd = mean(ACd,3);    acorrp = mean(ACp,3);
else
    mn    = M ;    mnd    = Md ;    mnp    = Mp ;
    sd    = SD;    sdd    = SDd;    sdp    = SDp;
    xcorr = XC;    xcorrd = XCd;    xcorrp = XCp;
    acorr = AC;    acorrd = ACd;    acorrp = ACp;
end

%% Figures

figure(1);
subplot(1,3,1); plot(pie); 
subplot(1,3,2); spy(Pstar); 
subplot(1,3,3); plot(a,Eq)

%% Default sets 

ds     = reshape(D,na,nex);
z2 =[min(z);z'; max(z); max(z)]; % Useful when filling
 
 for i=1:nex
hgt    = (sum(ds(:,i)));
yll(i) = (((amax-amin)/na)*hgt'+amin);
 end
yll  = yll';
yll2 = [amax; yll; yll(end); amax];
 
figure('color', 'w');
fill(yll2, z2, [0.9,0.9,0.9]);
xlim([amin amax]);
ylim([min(z) max(z)]);
set(gca,'Ytick',linspace(min(z),max(z),4),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Z','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

%% Tables

fprintf('\n\n\n      Sovereign Default Model \n') 
fprintf('\n                Parameters        \n') 
fprintf('\nDiscount factor                    %5.2f'   ,beta)
fprintf('\nRisk-free rate (pct)               %5.2f'   ,rf*100)
fprintf('\nRisk aversion                      %5.0f'   ,mu)
fprintf('\nProbability of reentry             %5.2f'   ,lambda)
fprintf('\nDefault penalty (pct)              %5.2f'   ,d0*100)
fprintf('\nBorrowing limit (pct)              %5.1f'   ,amin*100)
fprintf('\nLong run annual growth rate        %5.2f'   ,log(mn_g)*400)
fprintf('\nTransitory shock std dev (pct)     %5.2f'   ,sigma_z*100)
fprintf('\nTransitory shock autocorrelation   %5.2f'   ,rho_z)
fprintf('\n')
fprintf('\n           Ergodic Moments          \n') 
fprintf('\nNumber of iterations over q        %5.0f'   ,it)
fprintf('\nExpected Net Foreign Assets        %5.2f'   ,Eap)
fprintf('\nEq price of sovereign bond         %5.2f'   ,Eq(getindex(Eap,A)))
fprintf('\nErgodic default rate (pct)         %5.2f'   ,Dr*100)
fprintf('\n')
fprintf('\n')
fprintf('\n           Parameters and Data to match         \n') 
fprintf('\n Parameters         ') 
fprintf('\nBeta                        %5.3f'   ,beta)
fprintf('\nLambda                      %5.3f'   ,lambda)
fprintf('\nD0                          %5.3f'   ,d0)
fprintf('\nD1                          %5.3f'   ,d1)
fprintf('\n')
fprintf('\n                               Data          Model Unc.') 
fprintf('\n                         std.dev  acorr   std.dev  acorr') 
fprintf('\nCyclical Output          %5.2f   %5.2f    %5.2f   %5.2f'  , 4.18,0.96,sd(6),acorr(6))
fprintf('\nConsumption              %5.2f   %5.2f    %5.2f   %5.2f'  , 4.27,0.96,sdp(4),acorrp(4))
fprintf('\n')
fprintf('\n Data to match       ') 
fprintf('\n                            Data       Model        ') 
fprintf('\nDebt to GDP  (pct)          %5.2f      %5.2f' ,29.07    ,-mean(at(find((history==0)&(dt==0)))./yt((find((history==0)&(dt==0)))))*100)
fprintf('\nDefault rate (pct)          %5.3f      %5.3f' ,1.0      ,mn(5)*100)
fprintf('\nAverage EMBI (pct)          %5.3f      %5.3f' ,embi*100 ,(((1/mn(3))-1)-rf)*100)
fprintf('\nBond Price std. dev (pct)   %5.3f      %5.3f' ,1.84     ,sd(3))
fprintf('\n')
fprintf('\n           Simulated Moments          \n') 
fprintf('\n                                      Unconditional            If default           If no default')
fprintf('\n                                    mean  std.dev  corr    mean  std.dev  corr    mean  std.dev  corr')
fprintf('\nTotal Debt                         %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(1),sd(1),xcorr(1,2),mnd(1),sdd(1),xcorrd(1,2),mnp(1),sdp(1),xcorrp(1,2))
fprintf('\nBond Price                         %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(3),sd(3),xcorr(3,2),mnd(3),sdd(3),xcorrd(3,2),mnp(3),sdp(3),xcorrp(3,2))
fprintf('\nOutput                             %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(6),sd(6),xcorr(6,2),mnd(6),sdd(6),xcorrd(6,2),mnp(6),sdp(6),xcorrp(6,2))
fprintf('\nConsumption                        %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(4),sd(4),xcorr(4,2),mnd(4),sdd(4),xcorrd(4,2),mnp(4),sdp(4),xcorrp(4,2))
fprintf('\nCurrent Account                    %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(9),sd(9),xcorr(9,2),mnd(9),sdd(9),xcorrd(9,2),mnp(9),sdp(9),xcorrp(9,2))
fprintf('\nTrade Balance                      %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(10),sd(10),xcorr(10,2),mnd(10),sdd(10),xcorrd(10,2),mnp(10),sdp(10),xcorrp(10,2))
fprintf('\nEMBI                               %5.2f  %5.2f  %5.2f     %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f'   ,mn(11),sd(11),xcorr(11,2),mnd(11),sdd(11),xcorrd(11,2),mnp(11),sdp(11),xcorrp(11,2))
fprintf('\n')

%before_after_transition
