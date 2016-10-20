%% ENDOWMENT.M  Small Open Endowment Economy with Incomplete Markets
%
% Written by F. Hamann. Feel free to copy, change and distribute.
 fprintf('\nSmall open endowment economy with incomplete markets \n')

%% Parameters
 sigma	= 2;        % risk aversion
 beta	= 0.99;     % discount factor
 R      = 1.004;    % gross asset return rate (vs. R=1)

 if R*beta>=1; display('Set beta*R<1 for convergence'); break; end;

%% Markov chain for y
 [y,Py] = markovchain(2,0.5,0.8,0.5,0.5,0.8);

%% State-space S = YxB
 b = linspace(-15,10,500)';   
 
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

%% Simulation
 T      = 500;      
 s0     = findnearest(bmean,B);   
 spath  = simulmarkov(pstar,T,s0);

 ypath  = Y(spath);
 cpath  = ypath + R*B(spath)-b(x(spath));
 CApath = b(x(spath))-B(spath);

 sd.y = std(ypath);
 sd.c = std(cpath);

 figure(2)
 plot([ypath cpath])

 [sdev,corrcont,corr,acov] = samplemoms([ypath cpath CApath],1,3);

%% Model steady state statistics  
 fprintf('\nSteady-state Model Statistics \n ')
 fprintf('\nSample means ')
 fprintf('\n Earnings           %8.2f'  ,ymean)  
 fprintf('\n Consumption        %8.2f'  ,cmean) 
 fprintf('\n Net Assets         %8.2f'  ,bmean) 
 fprintf('\n Assets to income   %8.2f'  ,bmean/ymean) 
 fprintf('\nSample volatility') 
 fprintf('\n Consumption        %8.2f'  ,sd.c) 
 fprintf('\n Earnings           %8.2f\n',sd.y)

