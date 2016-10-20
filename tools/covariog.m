function c = covariog(X,maxtau);

% COVARIOG   Computes covariograms
% Usage:
%            c = covariog(X,maxtau)  
%
% calculates the covariogram for each column of the T x n matrix X where
%
%       c(tau,:)=cov( X(t,:), X(t-tau,:) ),  tau=-maxtau,...,maxtau
%
% Thus, the ith column of c is the covariogram for the ith column of X, i=1,...n.


[T,n]=size(X);
if T-maxtau<=0;
  str='The number of observations in the time series must be larger than tau';
  error(str);
end;
for i=1:n;
  for j=1:maxtau+1;
    vcov   = cov([X(j:T,i),X(1:T-j+1,i)]);
    c(j,i) = vcov(2,1);
  end;
end;
c = [c(maxtau+1:-1:2,:);c];
