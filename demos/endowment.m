%% ENDOWMENT.M  Small Open Endowment Economy with Incomplete Markets
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nSmall open endowment economy with incomplete markets \n')

%% Load data

 rgdp;                      % load Real GDP Colombia (1950-2014); FRED data

% FRED Graph Observations
% Federal Reserve Economic Data
% Link: https://fred.stlouisfed.org
% Help: https://fred.stlouisfed.org/help-faq
% Economic Research Division
% Federal Reserve Bank of St. Louis
% RGDPNACOA666NRUG 
% Real GDP at Constant National Prices for Colombia, Millions of 2011 USD, 
% Annual, Not Seasonally Adjusted

 lgdp       = log(RGDP)          ;
 lgdp_trend = hpfilter(lgdp,100) ;
 lgdp_cycle = lgdp-lgdp_trend    ;
 [rho,sige] = ols(lgdp_cycle(1:end-1),lgdp_cycle(2:end));
 [y,Py]     = rouwenhorst(9,0,rho,sige);
 y          = exp(y');

%% Model parameters
 sigma	= 2;           % risk aversion
 beta	= 0.9802;      % discount factor
 R      = 1.02;        % gross asset return rate (vs. R=1)

 if R*beta>=1; disp('Set beta*R<1 for convergence');  end;

%% State-space S = YxB
 b = linspace(-1,1,1000)';

 [Y,B] = gridmake(y,b);

 n = length(y)*length(b); 
 m = length(b);

%% Utility function and feasible consumption C>=0
 C = zeros(n,m);

 for i=1:m    
    C(:,i)=Y+R*B-b(i);  
 end

 C(C<=0) = NaN;
 u  = (C.^(1-sigma)-1)./(1-sigma);

%% Transition probability matrix (see Sargent and Ljundqvist)
 P = kron(speye(m,m),repmat(Py,m,1));

%% Bellman equation
 [v,x,pstar] = solvedp(u,P,beta,'policy');  clear P u C;

%% Steady State Distribution
 d = ergdist(pstar);
 
%% Summary statistics 
 c = Y+R*B-b(x);

 ymean = ergdist(Py)'*y;
 bmean = b(x)'*d;
 cmean = c'*d;

%% Plot some model properties
 plotdp(v,x,pstar,Y,B,y,b);

%% Model simulation
 T    = 500;      
 s0   = findnearest(bmean,B);   
 s_t  = simulmarkov(pstar,T,s0);

 y_t  = Y(s_t);
 c_t  = y_t + R*B(s_t)-b(x(s_t));
 CA_t = b(x(s_t))-B(s_t);

 sd.y = std(y_t);
 sd.c = std(c_t);

 figure; plotyy(1:T,y_t,1:T,c_t)
 
 figure; plot(CA_t)
  
 [sdev,corrcont,corr,acov] = samplemoms([y_t c_t CA_t],1,3);

%% Model steady state statistics  
 fprintf('\nSteady-state Model Statistics \n ')
 fprintf('\nSample means ')
 fprintf('\n Earnings           %8.2f'  ,ymean)  
 fprintf('\n Consumption        %8.2f'  ,cmean) 
 fprintf('\n Net Assets         %8.2f'  ,bmean) 
 fprintf('\n Assets to income   %8.2f'  ,bmean/ymean) 
 fprintf('\nSample volatility') 
 fprintf('\n Consumption        %8.3f'  ,sd.c) 
 fprintf('\n Earnings           %8.3f\n',sd.y)