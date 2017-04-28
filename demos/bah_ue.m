%% BAH_UE.M  Bewley-Aiyagari-Huggett vs. Uzawa-Epstein model
%            Follows Durdu, Mendoza and Terrones (2007) NBER-WP 13123
%
% Written by F. Hamann. Feel free to copy, change and distribute.

% % % Features
% 1) The Steady State values for Consumption and Foreing Assets are
% included
% 2) For the AR(1) discretization as a Markov chain it is included the
% Tauchen, G. & Hussey (1991) approach (as in the paper) coded by Martin
% Flodén. 
%   + Note: It is worth noting (and can be tested in the current model) that 
%           for high Autocorrelation values of the variable to be
%           approximated as a Markov chain, the Tauchen, G. & Hussey (1991) method and the
%           Flodén, M (2008) feature, don't perform very well in matching the
%           unconditional moments of the variable. See: Floden(2008). doi:10.1016/j.econlet.2007.09.040
%           Also, the Rouwenhorst (1995) method matches more accurately the unconditional moments than the
%           Tauchen-Hussey Method. See: Kopecky, K & Suen, R(2010) doi:10.1016/j.red.2010.02.002
% 3) It is constructed a grid for the BAH approach and another for the UE. Both grids differ 
%    in the imposed limits. For the BAH it is used the ad hoc borrowing
%    limit, for the UE it is used the natural borrowing limit.
% 4) It is modified the table with the results. The features in the table are: 
%    + The Markov Chain Approximation Statistics
%    + The Equivalent Precautionary Premium - EPP
%    + Coefficient of variation for each variable
%    + Output correlations
%    + First order autocorrelations for each variable.

% Note: It should be noted that the results in the present code may vary
% with respect to the article's results for two main reasons: one reason lies in the definition
% of the highest bound of the Foreign Assets' grid. The second reason, and
% maybe the most important one, is that the AR(1) discretization in the
% current code matches almost perfectly the autocorrelation of the
% variable, while in the article it is not as accurate.


 fprintf('\nDurdu, Mendoza and Terrones (2007) NBER-WP 13123')
 fprintf('\nSingle good SOE model calibrated for Mexico \n')

%% Parameters
 sigma	= 2;                         % risk aversion coefficient 2
%  sigma	= 5;                     % (Exercise 1) risk aversion coefficient 5

yss    = 1;                          % GDP trend mean            1

R      = 1.059;                      % gross asset return rate   1.059
phi    = -0.51;                      % ad-hoc borrowing limit    -0.51 

 
%% Data to match
 bss  = -0.44;                       % NFA to GDP ratio          -0.44
 css  = 0.692;                       % consumption to GDP        0.692
 abs  = yss+bss*(R-1)-css;           % absorption

 rhoy = 0.597;                       % GDP cycle autocorrelation 0.597
%  rhoy = 0.7;                       % (Exercise 2) GDP cycle autocorrelation 0.7

 sdvy = 0.03301;                     % GDP cycle std deviation   0.03301
%  sdvy = 0.05;                      % (Exercise 3) GDP cycle std deviation   0.05
%  sdvy = 0.025;                     % (Exercise 4) GDP cycle std deviation   0.025
  
%% Discount factor
 beta_bah = 0.94;                    % Bewley-Aiyagari-Huggett   0.94

 disc_ue  = log(R)/log(1+css);        % Uzawa-Epstein disc rate   
 beta_ue  = (1+css)^-disc_ue;         % Uzawa-Epstein disc factor
 
%% Uzawa-Epstein's Steady state
 css_ue = R^(1/disc_ue)-1;
 bss_ue = (css_ue-yss+abs)/(R-1);
 
%% Stationarity condition
 if R*beta_bah>=1; display('Set beta*R<1 for convergence'); end;

%% Markov chain for y

 sdve  = sqrt(sdvy^2*(1-rhoy^2));
 
% % % % With the Rouwenhorst(1995) method 
 p     = (1+rhoy)/2;
 q     = p;
 [y,Py] = markovchain(5,p,q,sdvy,yss);
 
% % With the Tauchen-Hussey (1991) method
%  omeg  = (1/2)+(rhoy/4);                  % Floden's Quadrature weighting value
%  sigmaZ = omeg*sdve + (1-omeg)*sdvy;	    % Floden's Quadrature 
%  [y,Py] = tauchenhussey(5,yss,rhoy,sdve,sigmaZ);

 spath  = simulmarkov(Py,20000,1); 
 
 ypath  = y(spath);
 
 [ymean,ysdv,ycorrcont,ycorr,acorry] = samplemoms(ypath,1,2);
 
%% State-space S = YxB

% % % BAH-approach's grid
%  b = linspace(phi,-phi,1500)';   


% % UE-approach's grid
 ndl = -min((y-abs)/(R-1));      % Natural Debt Limit for UE approach
%  The ad hoc borrowing constraint must be between the Natural Debt Limit
%  (ndl) and the Deterministic steady state (bss_ue). The limits are set to
%  place the ad hoc borrowing constraint (phi) in the middle of the grid.
 b = linspace(4*phi,-2*phi,1500)';   
 
 
 [Y,B] = gridmake(y,b);
 n = length(y)*length(b); 
 m = length(b);

%% Utility function and feasible consumption C>=0
 c = zeros(n,m);

 for i=1:m    
    c(:,i)=Y+R*B-b(i)-abs;  
 end

 c(c<=0) = NaN;
 u  = (c.^(1-sigma))./(1-sigma);
 
 beta = exp(-disc_ue*log(1+c));
 beta(beta>=1) = NaN;


%% Transition probability matrix
 P = kron(speye(m,m),repmat(Py,m,1));

%% Solve Bellman equation
%  [v,x,pstar] = solvedp(u,P,beta_bah,'policy');  % Bewley-Aiyagari-Huggett
 
[v,x,pstar] = solvedp(u,P,beta);                % Uzawa-Epstein
 clear P u C;
%% Steady State Distribution
 d = ergdist(pstar);
 
%% Summary statistics 
 c = Y+R*B-b(x)-abs;

 ymean = ergdist(Py)'*y;
 bmean = b(x)'*d;
 cmean = c'*d;
 sdc = sqrt(((c-cmean).^2)'*d);
 sdb  = sqrt(((b(x)-bmean).^2)'*d);

%% Plot some model properties
 plotdp(v,x,pstar,Y,B,y,b);

%% Simulation
 T      = 90000;      
 s0     = findnearest(bmean,B);   
 spath  = simulmarkov(pstar,T,s0);
 

 ypath  = Y(spath);
 cpath  = ypath + R*B(spath)-b(x(spath))-abs;
 CApath = b(x(spath))-B(spath);
 bpath  = b(x(spath));
 
 figure(2)
 plot([cpath ypath])
 
 [mn,sd,corr,acorr,auocorr] = samplemoms([ypath cpath CApath bpath],1,3);

  EPP    = (1+sigma)*(sd(2)/100/mn(2))^2/2*100;
  
  model = 1;                % write 'model = 1' to calculate precautionary
%       savings (ps) in the UE model using bss_ue. Otherwise the ps will be calculated using the BAH model phi.
  
  if model == 1
      phi = bss_ue;
  end
  
  ps     = bmean - phi;

  %% Markov chain approximation
  fprintf('\n   Markov Chain Approximation Statistics ')
 fprintf('\n                                  Approximation   Data %8.3f')
 fprintf('\n Mean                             %8.3f     %8.3f'  ,ymean,yss)  
 fprintf('\n Standard deviation (in percent)  %8.3f     %8.3f'  ,ysdv,sdvy*100)  
 fprintf('\n Autocorrelation                  %8.3f     %8.3f\n'  ,acorry(2),rhoy)  
 
 %% Model steady state statistics  
 fprintf('\n           Steady-state Model Statistics \n ')
 fprintf('\nPrecautionary savings (in percent)')
 fprintf('\n Prec. savings      %8.3f'    ,ps*100)  
 fprintf('\nEquivalent Precautionary Premium EPP ')
 fprintf('\n EPP                %8.3f'    ,EPP)  
 fprintf('\nSample means ')
 fprintf('\n Earnings           %8.3f'    ,mn(1))  
 fprintf('\n Consumption        %8.3f'    ,mn(2)) 
 fprintf('\n Net Assets         %8.3f'    ,bmean) 
 fprintf('\n Assets to income   %8.3f\n'  ,bmean/ymean) 
 fprintf('\nSample volatility (in percent)') 
 fprintf('\n Earnings           %8.3f'    ,sd(1))
 fprintf('\n Consumption        %8.3f'    ,sd(2)) 
 fprintf('\n Net Assets         %8.3f\n'  ,sdb*100) 
 fprintf('\nCoefficient of variation (in percent)') 
 fprintf('\n Earnings           %8.3f'    ,sd(1)/mn(1))
 fprintf('\n Net Assets         %8.3f'    ,sdb*100) 
 fprintf('\n Consumption        %8.3f\n'  ,sd(2)/mn(2)) 
 fprintf('\nOutput correlations') 
 fprintf('\n Earnings           %8.3f'    ,corr(1,4))
 fprintf('\n Consumption        %8.3f'    ,corr(2,4)) 
 fprintf('\n Current Account    %8.3f'    ,corr(3,4)) 
 fprintf('\n Net Assets         %8.3f\n'  ,corr(4,4)) 
 fprintf('\nAutocorrelations') 
 fprintf('\n Earnings           %8.3f'    ,auocorr(2,1))
 fprintf('\n Consumption        %8.3f'    ,auocorr(2,2)) 
 fprintf('\n Current Account    %8.3f'    ,auocorr(2,3)) 
 fprintf('\n Net Assets         %8.3f\n'  ,auocorr(2,4)) 
% toc
