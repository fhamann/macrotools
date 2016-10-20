function [Tran,s,probst,alambda,asigmay]=appmarkov2(const,lambda,sigma,m,N)

% APPMARKOV2 approximates a nonzero-mean AR(1) process with Markov chain
%
% Syntax: Let   y_t = const + lambda * y_(t-1) + u_t where
%         u_t: Gaussian white noise process with standard deviation sigma
% then
%        [Tran,s,probst,alambda,asigmay] = appmarkov(lambda,sigma,m,N)
% 
%    INPUTS:
%       lambda  = autocorrelation coefficient
%       sigma   = standard deviation of shocks
%       m       = "width" of discretized state space s
%       N       = number of possible states to approximate {y_t} (N=9 usually)
%
%    OUTPUTS:
%       Tran    = transition matrix of the Markov chain
%       s       = discretized state space of y_t
%       probst  = invariant distribution of Markov Chain
%       alambda = theoretical first order autoregression coeff. for Markov chain
%       asigma  = theoretical standard deviation for Markov chain Y_t
%
% NOTE: Inside the code the following options can be set:
%   tol      : convergence tolerance, default is sqrt(eps)
%   maxit    : maximum number of iterations, default is 100


maxit = 100;
tol   = sqrt(eps);

% Discretize the state space

stvy = sqrt(sigma^2/(1-lambda^2)); % standard deviation of y_t

% Equally space grid of points
% ymax=m*vary,ymin=-m*vary, ymax and ymin are two boundary points
ymax = m*stvy;                     % upper boundary of state space
ymin = -ymax;                      % lower boundary of state space
w = (ymax-ymin)/(N-1);             % length of interval 
s = ymin:w:ymax;                   % the discretized state space


% Calculate the transition matrix

for j=1:N;
   for k=2:N-1;
         Tran(j,k)= cdfnorm(s(k)-const-lambda*s(j)+w/2,0,sigma)...
         - cdfnorm(s(k)-const-lambda*s(j)-w/2,0,sigma);
   end
   Tran(j,1) = cdfnorm(s(1)-const-lambda*s(j)+w/2,0,sigma);
   Tran(j,N) = 1 - cdfnorm(s(N)-const-lambda*s(j)-w/2,0,sigma);   
end

if sum(Tran') ~= ones(1,N)
   str = find(abs(sum(Tran')-ones(1,N))>1e-6);  % find rows not adding up to one
   disp('Error in transition matrix');
   disp(['rows ',num2str(str),' does not sum to one']);
   disp('Normalizing the probabilities');
   d=sum(Tran');
   for i=1:length(str)
       Tran(str(i),:)=(1/d(str(i)))*Tran(str(i),:);
   end    
end


% Calculate the invariant distribution of Markov chain
Trans  = Tran';
probst = (1/N)*ones(N,1);            % initial distribution of states

for i=1:maxit
   probst1 = Trans*probst;
   if max(abs(probst1-probst))<tol, break, end;
   probst = probst1; 
end
   
   meanm = s*probst;                 % mean of invariant distribution
   varm = ((s-meanm).^2)*probst;     % variance of invariant distribution
    
   midaut1 = (s-meanm)'*(s-meanm);   % cross product of deviation from the
                                     % mean of y_t and y_t-1
                                   
   probmat = probst*ones(1,N);       % each column is invariant distribution   
      
   midaut2 = Tran.*probmat.*midaut1; % product of the first two terms is 
                                     % the joint distribution of (Y_t-1,Y_t)    
				     
   autcov1 = sum(sum(midaut2));      % first-order auto-covariance
  
   
% Calculate the asymptotic second moments of Markov chain
   alambda = autcov1/varm;            % theoretical lambda
   asigmay = sqrt(varm);