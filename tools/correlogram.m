function c = correlogram(X,p)

% CORRELOGRAM  Correlogram of a matrix of time series
%
% Usage:
%            c = correlogram(X,p)  
%
% Inputs:
%        X  : T by N matrix of N time series with T observations
%        p  : scalar, number of leads and lags for correlogram
%
% Output:
%        c  : 2*p+1 by N matrix of autocorrelations
%        
% Note: calculates correlogram for each column of the TxN matrix X where
%
%       c(s,:)=corr(X(t+s,:),X(t-s),:) ),  s = -p,...,p
%
% The ith column of c is the correlogram for the ith column of X, i=1,...N
% Also note that position s=(p+1), c(s,:)=corr(X(t+s,:),X(t-s,:)) = 1.

[T,n]=size(X);
if T-p<=0
  str='Length of time series must be larger than the number of input lags';
  error(str);
end
for i=1:n
  for j=1:p+1
    vcov   = nancov([X(j:T,i),X(1:T-j+1,i)]);
    c(j,i) = vcov(2,1)/((sqrt(vcov(2,2))*sqrt(vcov(1,1))));
  end
end
c = [c(p+1:-1:2,:);c];