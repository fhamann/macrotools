function c = correlogram(X,maxtau)

% CORRELOGRAM   Correlograms
%
% Usage:
%            c = correlogram(X,maxtau)  
%
% Note: calculates correlogram for each column of the Txn matrix X where
%
%       c(tau,:)=corr( X(t,:), X(t-tau,:) ),  tau=-maxtau,...,maxtau
%
% The ith column of c is the correlogram for the ith column of X, i=1,...n


[T,n]=size(X);
if T-maxtau<=0;
  str='Number of observations in time series must be larger than tau';
  error(str);
end;
for i=1:n;
  for j=1:maxtau+1;
    vcov   = cov([X(j:T,i),X(1:T-j+1,i)]);
    c(j,i) = vcov(2,1)/((vcov(2,2)*vcov(1,1))^0.5);
  end;
end;
c = [c(maxtau+1:-1:2,:);c];
