% OIL_EG 
%
% Stochastic oil price and endowment & sovereign default
% Sovereign owns and operates oil extraction company decides optimal oil 
% extraction, borrowing and default policies.
% Default decision as in Eaton and Gersovitz with exogenous reentry

close all; clear all;

%% Model parameters  
  beta   = .95;             % Discount factor                           .95
  mu     = 2;               % KEY: Risk aversion                      {2-5}
  kappa  = 10;              % KEY: Mg cost of extraction             {1-15}
  omax   = 36;              % max oil reserves                           36
  omin   = 31;              % min oil reserves                           31
  K      = 1;               % discoveries                                 1
  Ep     = 1;               % expected relative oil price                 1
  Ey     = 1;               % expected endowment of consumption good      1
  sigmay = .05;             % KEY: sdv shocks                           .05
  sigmap = .3;              % KEY: sdv shocks                            .3
  rho    = 0.5;             % auto-correlation of shocks                 .5
  rhopy  = 0;               % correlation oil price - endowment           0
  Rfree  = 1.02;            % risk free rate                           1.02
  delta  = 0.0;             % output loss in autarky                      0
  lambda = 0.1;             % prob of redemption                        0.1
  bmin   = -.5*Ey;          % borrow. limit (% of non-oil output)       -.5
  bmax   = 0;               % lending to the ROW                          0
 
%% Construct action and state space
% Build the space of states: the sets "O" and "B"

  o      = (omin:.3:omax)'          ;  no = length(o) ;
  b      = linspace(bmin,bmax,20)'  ;  nb = length(b) ;
 
% Build the space of ACTIONS: OxB. That is, the next periods value of the 
% endogenous state variables. 
% Note: "op" denotes "o prime" and "bp" denotes "b prime"

  [op,bp] =  gridmake(o,b)           ;  m  = length(op) ;  
      
% Build the space of STATES: OxBxPxY

  y = Ey*exp([-sigmay;sigmay])       ;
  p = Ep*exp([-sigmap;sigmap])       ;  nz = length([p; y]);
  
  [O,B,P,Y] = gridmake(o,b,p,y)      ;  n  = length(O);
  
%% Set probability transition matrix
% See Mendoza (1991) for the setup of the 2-state 2-point Markov Chain.
% Two exogenous state variables: p and y
% Two possible values: low and high

  Ppy  = ptm(rho,rhopy);
 
%% Declare consumption and profits to store feasible values

 cpay   = zeros(n,m)          ;       cdef   = zeros(n,m);
 piepay = zeros(n,m)          ;       piedef = zeros(n,m);

 bzero  = find(b==0)          ;       % state position where debt is zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outer-loop: find "q" that satisfies the capital market equilibrium     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxit = 10                    ;       % maximum number of iterations
tol   = 10e-6                 ;       % tolerance

% Step 1: Guess an initial q

q    = (1/(Rfree))*ones(n,m)  ;       % guessed bond price (risk-free)
R    = 1./q                   ;       % guessed gross interest rate 

for it = 1:maxit
  qold = q;    
  Rold = R;
  % 1. Construct the reward function (with and w/o default), for a given q 
  
   for i=1:n
    for j=1:m
     if O(i)-op(j)+K>=0 && O(i)-op(j)+K<=O(i)
        % Profits under repayment state with LINEAR (in x) cost function
        piepay(i,j) = P(i)*(O(i)-op(j)+K)-kappa*((O(i)-op(j)+K)^2/(O(i)));
        % Profits under default state with LINEAR (in x) cost function
        piedef(i,j) = P(i)*(O(i)-op(j)+K)-kappa*((O(i)-op(j)+K)^2/(O(i)));             
     else
         piepay(i,j) = -inf;
         piedef(i,j) = -inf;       
     end  
        % Consumption under repayment state      
         cpay(i,j)   = Y(i) + piepay(i,j) + B(i) - q(i,j)*bp(j);
        % Consumption under default state          
         cdef(i,j)   = (1-delta)*(Y(i) + piedef(i,j));
    end
   end
    
   upay = (cpay.^(1-mu)-1)./(1-mu); upay(cpay<=0)=-Inf; 
   udef = (cdef.^(1-mu)-1)./(1-mu); udef(cdef<=0)=-Inf;
   
  % 2. Solve the optimal stopping problem, for a given q
  
  [vd,xd,vr,xr,default] = solvedpos(upay,udef,Ppy,beta,lambda,bzero);

  % 3. Compute the expected value of default, for a given q

  Edef = Ppy*(reshape(default,m,nz)');
  ind  = reshape(default,m,nz)';
  ind  = sum(ind)==nz;                 % find b' where default is Prob=1 

  for j=1:m;
      if ind(j)==1;
        Edef(:,j)=1;
      end;
  end;

  Edef = kron(Edef,ones(m,1)); % add b(t) to rows to make same size as q

  % 4. Update the price "q" (i.e. find a new "q")

  q = (1/Rfree)*(1-Edef);
  q = max(q,0);
  R = 1./q;

  % 5. Stop if the old "q" is similar to new "q"

  if abs(qold-q)<tol,break,end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Outer loop: already found the "q"                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optimal policy & transition prob matrix, for the equilibrium "q"
% [vd,xd,vr,xr,default,Pr] = solvedpos(upay,udef,Ppy,beta,lambda,bzero,2,1);

x       = (1-default).*xr+default.*bzero;                % policy function
% plotyy(1:n,policy,1:n,[vd vr]);

posGrid = gridmake((1:m)',(1:length(Ppy))');             % states positions
Pr      = otpmos(x,Ppy,lambda,bzero,vd,vr,posGrid,2,1);  % OTPM

%% Ergodic distribution of oil reserves and debt, for the equilibrium "q"

pie  = ergdist(Pr);    

%% Makeplots

F=[P, Y];
[r1, ~]=find(F(:,1)==p(1) & F(:,2)==y(1)); % find Plow and Ylow positions
[r2, ~]=find(F(:,1)==p(1) & F(:,2)==y(2)); % find Plow and Yhigh positions
[r3, ~]=find(F(:,1)==p(2) & F(:,2)==y(1)); % find Phigh and Ylow positions
[r4, ~]=find(F(:,1)==p(2) & F(:,2)==y(2)); % find Phigh and Yhigh positions

vd_ll=reshape(vd(r1), [no nb]); vd_lh=reshape(vd(r2), [no nb]);
vd_hl=reshape(vd(r3), [no nb]); vd_hh=reshape(vd(r4), [no nb]);
vr_ll=reshape(vr(r1), [no nb]); vr_lh=reshape(vr(r2), [no nb]);
vr_hl=reshape(vr(r3), [no nb]); vr_hh=reshape(vr(r4), [no nb]);


pie_ll=reshape(pie(r1), [no nb]); pie_lh=reshape(pie(r2), [no nb]);
pie_hl=reshape(pie(r3), [no nb]); pie_hh=reshape(pie(r4), [no nb]);

% Figure 1: P low and Y low
figure('color', 'w');subplot(2,2,1);
surf(b,o, vd_ll);
hold on
mesh(b,o,vr_ll);
legend('$V_{default}$','$V_{repay}$','Location','northeast')
title(['Value under default and repayment, with $P_{LOW}, Y_{LOW}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)
hold off

% Figure 2: P low and Y high
subplot(2,2,2)
surf(b,o, vd_lh);
hold on
mesh(b,o,vr_lh);
legend('$V_{default}$','$V_{repay}$','Location','northeast')
title(['Value under default and repayment, with $P_{LOW}, Y_{HIGH}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)
hold off

% Figure 3: P high and Y low

subplot(2,2,3);
surf(b,o, vd_hl);
hold on
mesh(b,o,vr_hl);
legend('$V_{default}$','$V_{repay}$','Location','northeast')
title(['Value under default and repayment, with $P_{HIGH}, Y_{LOW}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)
hold off

% Figure 4: P high and Y high

subplot(2,2,4); surf(b,o, vd_hh);
hold on
mesh(b,o,vr_hh);
legend('$V_{default}$','$V_{repay}$','Location','northeast')
title(['Value under default and repayment, with $P_{HIGH}, Y_{HIGH}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)
hold off

%Default and Repayment sets

ds_ll=(vd_ll==vr_ll); % Decision set with P low and Y low
ds_lh=(vd_lh==vr_lh); % Decision set with P low and Y high
ds_hl=(vd_hl==vr_hl); % Decision set with P high and Y low
ds_hh=(vd_hh==vr_hh); % Decision set with P high and Y high

b2=[bmin;b; bmax];    % Useful when filling

%Figure 1: Decision set with P low and Y low

for i=1:nb
hgt=(sum(ds_ll(:,i))-0);
yll(i)=(((omax-omin)/no)*hgt'+omin); 
end 
yll=yll';
yll2=[omax; yll; omax];

figure('color', 'w');subplot(2,2,1);
fill(b2, yll2, [0.9,0.9,0.9]);
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
ylim([omin omax]);
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment, with $P_{LOW}, Y_{LOW}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

%Figure 2: Decision set with P low and Y high

for i=1:nb
hgt=(sum(ds_lh(:,i))-0);
ylh(i)=(((omax-omin)/no)*hgt'+omin); 
end 
ylh=ylh';
ylh2=[omax; ylh; omax];

subplot(2,2,2); fill(b2, ylh2, [0.9,0.9,0.9])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
ylim([omin omax]);
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment, with $P_{LOW}, Y_{HIGH}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

%Figure 3: Decision set with P high and Y low

for i=1:nb
hgt=(sum(ds_hl(:,i))-0);
yhl(i)=(((omax-omin)/no)*hgt'+omin); 
end 
yhl=yhl';
yhl2=[omax; yhl; omax];

subplot(2,2,3); fill(b2, yhl2, [0.9,0.9,0.9])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
ylim([omin omax]);
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment, with $P_{HIGH}, Y_{LOW}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

%Figure 4: Decision set with P high and Y high

for i=1:nb
hgt=(sum(ds_hh(:,i))-0);
yhh(i)=(((omax-omin)/no)*hgt'+omin); 
end 
yhh=yhh';
yhh2=[omax; yhh; omax];

subplot(2,2,4); fill(b2, yhh2, [0.9,0.9,0.9])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
ylim([omin omax]);
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment, with $P_{HIGH}, Y_{HIGH}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')



% Figure 1: Ergodic distribution conditioned on P low and Y low
figure('color', 'w');subplot(2,2,1);
surf(b,o, pie_ll);
title(['Ergodic distribution conditioned on $P_{LOW}, Y_{LOW}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)

% Figure 2: Ergodic distribution conditioned on P low and Y high
subplot(2,2,2)
surf(b,o, pie_lh);
title(['Ergodic distribution conditioned on $P_{LOW}, Y_{HIGH}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)

% Figure 3: Ergodic distribution conditioned on P high and Y low

subplot(2,2,3);
surf(b,o, pie_hl);
title(['Ergodic distribution conditioned on $P_{HIGH}, Y_{LOW}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)

% Figure 4: Ergodic distribution conditioned on P high and Y high

subplot(2,2,4); 
surf(b,o, pie_hh);
title(['Ergodic distribution conditioned on $P_{HIGH}, Y_{HIGH}$.'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Oil reserves $o_{t}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
xlim([bmin bmax])
ylim([omin omax])
set(gca,'Xtick',linspace(bmin,bmax,6),'fontsize', 11)
set(gca,'Ytick',linspace(omin,omax,6),'fontsize', 11)
