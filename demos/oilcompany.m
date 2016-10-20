%% OILCOMPANY  Solves optimal oil extraction problem by private firm 

%% Model parameters  

Ep     = 1;          % expected oil price (normalized to 1)
rho    = 0;          % auto-correlation of oil price shocks (0 is iid)                     
sigma  = 0.2;        % sdv price shocks (0.15->25% std dev oil price)     
q      = .9645;      % disc rate (0.95 -> about 5% real interest rate)
kappa  = 2.5;        % mg cost extraction parameter 2.5
smax   = 1;          % possible oil reserves (normalized to 1)
d      = 0.05;       % fixed costless discoveries
fc     = 0*d;        % fixed cost

%% Construct state space PxS

s      = (0.0001:d/20:smax)';              ns = length(s);  

[p,Pp] = rouwenhorst(9,Ep,rho,sigma);      np = length(p);    
[P,S]  = gridmake(p',s);                   n  = length(S);

%% Profit function and feasible extraction x<=o
x   = zeros(n,ns);
pie = zeros(n,ns);
 
for i=1:ns*np
    for j=1:ns
        x(i,j) = (S(i)-s(j)+d);
        if x(i,j)>=0 && x(i,j)<=S(i)+d           
          pie(i,j) = P(i)*x(i,j)-kappa*((x(i,j))^2/(S(i)))-fc;
        else
          pie(i,j) = -inf;
        end
    end
end
  
%% Transition probability matrix (see Sargent and Ljundqvist)  
  
Pr = kron(speye(ns,ns),repmat(Pp,ns,1));

%% Solve Bellman equation for oil-firm

[v,x,pstar] = solvedp(pie,Pr,q,'policy');  clear Pr u c;
  
X = S + d - s(x);   % optimal extraction

%% Steady State Density, f

f = ergdist(pstar);

%% Compute optimal steady-state oil reserves desired by oil company
  
Eo  = f'*s(x); 
Ex  = f'*X;
  
profit = f'*(P.*X-kappa/2*(X.^2)./(1+S)-fc); 
value  = f'*(P.*S);
    
fprintf('Steady State Means\n') 
fprintf('   Stock           = %5.2f\n'  ,Eo)
fprintf('   Value of Stock  = %5.2f\n'  ,value)  
fprintf('   Extraction      = %5.2f\n'  ,Ex/Eo)
fprintf('   Yrs of Reserves = %5.2f\n'  ,Eo/Ex)
fprintf('   Profits         = %5.2f\n'  ,profit)  

%% Some plots 
% Reshaping
vr = reshape(v,np,ns)';
Xr = reshape(X,np,ns)';
fr = reshape(f,np,ns)';
Fr = cumsum(fr);
  
% Plot optimal policy
figure(1); 
  h=plot(s,Xr); % set(h,'FaceColor',[.75 .75 .75])
  axis([0 max(s) -inf inf]);
  title('Optimal Extraction Policy');
  xlabel('Oil Reserves'); ylabel('Oil Extraction');
  if np==2; legend('Low oil prices','High oil prices'); end

% Plot optimal value function
figure(2); plot(s,vr); 
  title('Optimal Value Function');
  xlabel('Oil Reserves'); ylabel('Value');
  if np==2; legend('Low oil prices','High oil prices'); end

% Compute steady state distribution of oil reserves
figure(3); 
  h=plot(fr); 
  title('Steady State Density');
  xlabel('Oil Reserves'); ylabel('Probability'); 
  if np==2; legend('Low oil prices','High oil prices'); end
  
% Compute steady state distribution of oil reserves
figure(4); 
  h=plot(Fr); 
  title('Steady State Distribution');
  xlabel('Oil Reserves'); ylabel('Probability'); 
  if np==2; legend('Low oil prices','High oil prices'); end
   