function [smean,sdev,corr,acorr] = samplemoms(x,m,k) 

% SAMPLEMOMS Computes time series moments
%
% Usage: 
%    [smean,sdev,corrcont,acorr] = samplemoms(x,m,k)
%
% Inputs:
%  x      : TxN matrix of N time series column-vectors sized T by 1 each
%  m      : position of pivot time series vector for cross-correlogram
%  k      : number of lags for cross-correlogram
% Ouputs:
%  smean : Nx1 vector of mean of N time series vectors
%  sdev  : Nx1 vector of standard deviations of N time series vectors
%  corr  : Nxk matrix of cross-colegrograms with respect to m for k lags
%  acorr : correlogram of the time series vectors contained in matrix x
%  


[nobs,nseries]=size(x);
if nargin<3; k = 3; end         % default number of k lags and k leads  

smean  = mean(x)';
sdev   = 100*std(x)';
corr   = ones(nseries,2*k+1);

for i=-k:k,
   if i<=0
      aux = corrcoef([x(1-i:nobs,m) x(1:nobs+i,1:nseries)]);
      corr(:,k+i+1) = aux(2:nseries+1,1);
  elseif i>0
      aux = corrcoef([x(1:nobs-i,m) x(1+i:nobs,1:nseries)]);
      corr(:,k+i+1) = aux(2:nseries+1,1);
   end
end

[acorr] = correlogram(x,k);