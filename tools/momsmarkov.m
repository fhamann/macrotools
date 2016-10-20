function [mmean,msdev] = momsmarkov(X,d)

% MOMSMARKOV Moments of controlled Markov chain X with ergodic distribution d
%
% Usage:  [mmean,msdev] = momsmarkov(X,d)
%
% Input:
%  X      : n by N matrix of N column-vectors sized n by 1 each, where
%           n is the number of discrete states and N is the number of
%           variables
%  d      : stationary density of controlled Markov process
%
% Ouput:
%  mmean  : Nx1 vector of mean of N time series vectors
%  msdev  : Nx1 vector of standard deviations of N time series vectors

[n,N] = size(X);
mmean = X'*d;
msdev = zeros(N,1);

for i=1:N
    msdev(i)=sqrt(((X(:,i)-mmean(i)).^2)'*d);
end