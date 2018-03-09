function ta = correlogram2(y,p)

% CORRELOGRAM2 - Compute Autocorrelations through p Lags
%  
%
% Inputs:
%       y : series to compute acf for, nx1 column vector
%       p : total number of lags, 1x1 integer
%
% Output:
%       c : px1 vector containing autocorrelations
%           (First lag computed is lag 1. Lag 0 not computed)
%

[n1, n2] = size(y) ;
if n2 ~=1
    error('Input y must be an nx1 vector. For matrix, use "correlogram".')
end

[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
error('Length of time series must be larger than the number of input lags')
end


ta = zeros(p,1) ;
N = max(size(y)) ;
ybar = mean(y); 

% Collect ACFs at each lag i
for i = 1:p
   ta(i) = acf_k(y,i) ; 
end


function ta2 = acf_k(y,k)
 
N    = max(size(y)) ;
ybar = mean(y); 
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) ;
end

% Denominator, unscaled variance
yvar = (y-ybar)'*(y-ybar) ;

ta2 = sum(cross_sum) / yvar ;
