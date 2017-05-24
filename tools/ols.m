function [beta,sdu,se,u,r2,rb2]=ols(X,y)

%OLS   Ordinary least-square estimation
%
% Usage:      [beta,sdu,se,u,r2,rb2]=ols(X,y)
%
%      calculates the OLS estimates for:
% 
%           y = X*beta+u,   E[u]=0, E[u*u']=sigma^2 * I
%     
%      where y is a Txn matrix of T observations on the n dependent 
%      variables, X is a Txk matrix of T observations on k independent
%      variables, and u is a Txn matrix of residuals.  If requested, 
%      OLS also returns the std deviation errors of u, sdu, the estimated 
%      residuals, u, the R^2 statistic, the Rbar^2 statistic and the
%      standard errors of residuals, se.


beta=(X'*X)\(X'*y);
if nargout>1;
  u     = y-X*beta;
  [T,k] = size(X);
  df    = T-k;
  sdu   = sqrt(diag(u'*u)/df);
  se    = sqrt(diag(inv(X'*X)/df))*sqrt(diag(u'*u))';
end;
if nargout>3;
  p   = X*((X'*X)\X');
  l   = eye(T)-ones(T)/T;
  r2  = ((diag(y'*X*beta)-(sum(y)').^2/T).^2)./(diag(y'*l*y).*diag(y'*p*l*p*y));
  rb2 = (T*r2-k)/df;
end;
