function [X,P] = markovchain(N,p,q,e,m,a)

% MARKOVCHAIN Markov process and transition matrix 
% 
% Purpose: Follows Rouwenhorst to find a Markov chain whose sample paths
%          approximate a family of random walk processes with mean m, std
%          dev e and serial correlation p+q-1. Depending on p and q, the
%          process could be AR(1) or ARCH.
%          See: Cooley's RBC book, chapter 10 by K. G. Rouwenhorst
%
% Usage:  [X,P] = markovchain(N,p,q,e,m,a)
%
% Input:   N   scalar, number of nodes for X, usually odd
%          p   scalar, governs cond. prob of moving up after low state
%          q   scalar, governs cond. prob of moving down after low state
%          e   scalar, size of epsilons
%          Optional: 
%          m   scalar, value of mean for random variable
%          a   scalar, value of downside asymmetry
%
% Output:  X   N*1 vector, nodes for X
%          P   N*N matrix, transition probabilities
%
% NOTE: For an evenly spaced grid of X:
%        - if 0<p=q<1 then X follows a standard AR(1) process.
%        - if p=q=0.5 then X has no persistence.
%        - if p<q bad times are longer than good ones (ARCH process).
%
%       Also, note that in the last case, altering p and q does not 
%       preserve the mean of X (i.e. is not a mean-preserving shift).
%
% Code by Franz Hamann, Banco de la Republica, December 2013
% Addional comments welcome


 P  = [p 1-p; 1-q q];

 for n = 3:N
    P = p*[P zeros(n-1,1); zeros(1,n)]      + ...
         (1-p)*[zeros(n-1,1) P; zeros(1,n)] + ...
         (1-q)*[zeros(1,n); P zeros(n-1,1)] + ...
         q*[zeros(1,n); zeros(n-1,1) P]     ;
    P(2:end-1,:) = P(2:end-1,:)/2;
 end

 z = sqrt(N-1)*e;
 X = linspace(-z,z,N)';
 
 if nargin>4;  X    = m+X   ; end   % adds mean, if requested
 if nargin>5;  X(1) = m-a*z ; end   % asymmetry, if requested  