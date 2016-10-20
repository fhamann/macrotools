function [BigPI,S,Q] = garch(b0,b1,b2,theta,lambda,p_0,I_p,q_1,I_q,m,n)

% GARCH approximates a GARCH(1,1) process with Markov chain
%
% Model:    x_t = p_0 - 1/2*h_t +sqrt(h_t)*u_t 
%           h_t = b0 + b1*h_(t-1) + b2*h_(t-1)*(u_(t-1)-theta-lambda)^2
%       where
%             u(t+1) = e(t+1) + lambda 
%       is normal standard zero-mean noise process.
%
% Syntax: [X,PI] = garch(b0,b1,b2,theta,lambda,p_0,m,n,T,tao)
% 
%    INPUTS:
%       b0     : scalar, parameter
%       b1     : scalar, parameter
%       b2     : scalar, parameter
%       theta  : scalar, leverage parameter
%       lambda : scalar, parameter
%       m      : scalar, number of discrete values for x_t
%       n      : scalar, number of discrete values for h_t
%       p_0    : scalar, constant of x_t equation
%       I_p    : scalar, width of x_t centered around p_0
%       q_1    : scalar, constant of h_t equation
%       I_q    : scalar, width of h_t centered around q_1
%
%    OUTPUTS:
%       BigPI  : transition matrix of the Markov chain
%       S      : discretized state space of x_t
%       Q      : discretized state space of h_t

% Discretiza GARCH-M como en: Duan & Simonato (2001) "American option 
% pricing under GARCH by a Markov Chain approx" JEDC, pp 1689-1718.

if b1+b2*(1+(theta+lambda)^2) >= 1
  'beta_0 beta_1 beta_2 theta lambda'  
  [b0 b1 b2 theta lambda]
  error('Nonstationary variance: See Duan & Simonato (2001), pp 1708.');
end

% Discretize the state space
h_star = b0/(1-b1-b2*(1+(theta+lambda)^2)); % definida en pag 1693 entre eq (6) y eq (7) ;

pbar = linspace(p_0-I_p,p_0+I_p,m);
qbar = linspace(q_1-I_q,q_1+I_q,n); % 
c = zeros(m+1); c(1)=-Inf; c(m+1)=Inf;
for i=2:m
    c(i)= (pbar(i)+pbar(i-1))/2; % pag 1698 abajo de eq (20);
end
d=zeros(n+1); d(1)=-Inf;d(n+1)=Inf;
for j=2:n
    d(j)=(qbar(j)+qbar(j-1))/2; % pag 1698 abajo de eq (21);
end
PI=zeros(m,n,m,n);
for k=1:m
 for l=1:n
  for i=1:m
   for j=1:n
    Phi(j,k,i) = log(b0+b1*exp(qbar(j))+b2*(pbar(k)-pbar(i)+...
                 0.5*(exp(qbar(j))-h_star)-(theta+lambda)*exp(qbar(j)/2))^2);               % eq (22) ;
    Lk(i,j)    = (c(k)  -pbar(i)+0.5*(exp(qbar(j))-h_star))/sqrt(exp(qbar(j)));             % eq (25) ;
    Lk1(i,j)   = (c(k+1)-pbar(i)+0.5*(exp(qbar(j))-h_star))/sqrt(exp(qbar(j)));             % eq (25) ;
    PI(i,j,k,l) = (normcdf(Lk1(i,j),0,1) - normcdf(Lk(i,j),0,1))*(d(l)<Phi(j,k,i)<=d(l+1)); % eq (23) ;
   end
  end
 end
end
% PI debe ser la matriz en pag 1695 eq (13) ;
BigPI = zeros(n*m,n*m); 
fila=1;
for j=1:n
    for i=1:m
        BigPI(fila,:)= reshape(squeeze(PI(i,j,:,:)),1,m*n);
        fila=fila+1;
    end
end

S = kron(ones(n,1),pbar');
Q = kron(qbar',ones(m,1));