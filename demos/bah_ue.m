%% BAH_UE.M  Bewley-Aiyagari-Huggett vs. Uzawa-Epstein model
%            Follows Durdu, Mendoza and Terrones (2007) NBER-WP 13123
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nDurdu, Mendoza and Terrones (2007) NBER-WP 13123')
 fprintf('\nSingle good SOE model calibrated for Mexico \n')

%% Parameters
 sigma	= 2;                         % risk aversion coefficient 2
 yss    = 1;                         % GDP trend mean            1
 R      = 1.059;                     % gross asset return rate   1.059
 phi    = -0.51;                     % ad-hoc borrowing limit    -0.51
 
%% Data to match
 bss  = -0.44;                       % NFA to GDP ratio          -0.44
 css  = 0.692;                       % consumption to GDP        0.692
 abs  = yss+bss*(R-1)-css;           % absorption
 rhoy = 0.597;                       % GDP cycle autocorrelation 0.597
 sdvy = 0.03301;                     % GDP cycle std deviation   0.03301
  
%% Discount factor
 beta_bah = 0.94;                    % Bewley-Aiyagari-Huggett   0.94
 disc_ue  = log(R)/log(1+css)        % Uzawa-Epstein disc rate   
 beta_ue  = (1+css)^-disc_ue         % Uzawa-Epstein disc factor

%% Stationarity condition
 if R*beta_bah>=1; display('Set beta*R<1 for convergence'); break; end;

%% Markov chain for y
 sdve  = sqrt(sdvy^2*(1-rhoy^2))
 p     = (1+rhoy)/2;
 q     = p;

 [y,Py] = markovchain(2,p,q,sdvy,yss);

 spath  = simulmarkov(Py,50000,1); ypath  = y(spath);
 
 [ymean,ysdv,ycorrcont] =samplemoms(ypath,1,2)
 
%% State-space S = YxB
 b = linspace(phi,0,1000)';   
 
 [Y,B] = gridmake(y,b);

 n = length(y)*length(b); 
 m = length(b);

%% Utility function and feasible consumption C>=0
 c = zeros(n,m);

 for i=1:m    
    c(:,i)=Y+R*B-b(i)-abs;  
 end

 c(c<=0) = NaN;
 u  = (c.^(1-sigma)-1)./(1-sigma);
 
 beta = exp(-disc_ue*log(1+c));
 beta(beta>=1) = NaN;


%% Transition probability matrix
 P = kron(speye(m,m),repmat(Py,m,1));

%% Solve Bellman equation
 [v,x,pstar] = solvedp(u,P,beta_bah,'policy');  % Bewley-Aiyagari-Huggett
% [v,x,pstar] = solvedp(u,P,beta);                % Uzawa-Epstein
 clear P u C;
%% Steady State Distribution
 d = ergdist(pstar);
 
%% Summary statistics 
 c = Y+R*B-b(x)-abs;

 ymean = ergdist(Py)'*y;
 bmean = b(x)'*d;
 cmean = c'*d;

%% Plot some model properties
 plotdp(v,x,pstar,Y,B,y,b);

%% Simulation
 T      = 500;      
 s0     = findnearest(bmean,B);   
 spath  = simulmarkov(pstar,T,s0);

 ypath  = Y(spath);
 cpath  = ypath + R*B(spath)-b(x(spath))-abs;
 CApath = b(x(spath))-B(spath);

 figure(2)
 plot([ypath cpath])

 [mn,sd,corr,acorr] = samplemoms([ypath cpath CApath],1,3)

%% Model steady state statistics  
 fprintf('\nSteady-state Model Statistics \n ')
 fprintf('\nSample means ')
 fprintf('\n Earnings           %8.3f'  ,mn(1))  
 fprintf('\n Consumption        %8.3f'  ,mn(2)) 
 fprintf('\n Net Assets         %8.3f'  ,bmean) 
 fprintf('\n Assets to income   %8.3f'  ,bmean/ymean) 
 fprintf('\nSample volatility') 
 fprintf('\n Consumption        %8.3f'  ,sd(2)) 
 fprintf('\n Earnings           %8.3f\n',sd(1))

