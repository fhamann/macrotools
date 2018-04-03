function [smean,sdev,corr,acorr] = samplemoms(x,m,k) 

% SAMPLEMOMS Time series moments of matrix of data
%
% Usage: 
%          [smean,sdev,corrcont,acorr] = samplemoms(x,m,k)
%
% Inputs:
%  x      : TxN matrix of N time series column-vectors sized T by 1 each
%  m      : position of pivot time series vector for cross-correlogram
%  k      : number of lags for cross-correlogram
%
% Ouputs:
%  smean : Nx1 vector, each element is mean of each N time series
%  sdev  : Nx1 vector, each element is std deviation of each N time series
%  corr  : Nxk matrix of cross-colegrograms with respect to m for k lags
%  acorr : correlogram of the time series vectors contained in matrix x 

[nobs,nseries] = size(x);

if nargin<3; k = 3; end         % default number of k lags and k leads  

smean  = nanmean(x)';
sdev   = nanstd(x)';
corr   = ones(nseries,2*k+1);

for i = 1:nseries
auxcorr   = xcorrelogram(x(:,m),x(:,i),k);
corr(i,:) = auxcorr';
end

for j = 1: nseries  
        aux1(:,j) = correlogram2(x(:,j),k);
end

acorr = [flip(aux1,1);ones(1,nseries);aux1];