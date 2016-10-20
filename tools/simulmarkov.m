function S = simulmarkov(P,T,S0);

% SIMULMARKOV  simulation of Markov chain with transition matrix P
%
% Purpose:   Follows Rouwenhorst to find a Markov chain whose sample paths
%            approximate a family of random walk processes with mean m, std
%            dev e and serial correlation p+q-1. Depending on p and q, the
%             process could be AR(1) or ARCH.
%            See: Cooley's RBC book, chapter 10 by K. G. Rouwenhorst
%
% Usage:     S = simulmarkov(P,T,S0)
%
% Input:     Q   NxN matrix, transition probabilities
%            T   scalar, number of periods to simulate
%            S0  scalar, initial state 
%
% Output:    S   Tx1 vector, simulated state path
%
%  S = simulmarkov(P,T) returns a T x 1 vector S, which is a realization
%  of a stationary Markov chain that has transition probabilities
%  in P.  The elements of S are integers from 1 to n where n is
%  the dimension of P.
%
%  S = simulmarkov(P,T,S0) uses the initial state S0 for S(1).

%  Code by    : Ellen McGrattan, 9-7-92 [Revised, ERM, 5-23-94]
%  Adapted by : Franz A Hamann , 3-5-02

n = length(P);

if nargin<3;
  S(1) = ceil(rand*n);
else
  S(1) = S0;
end;

F = cumsum(P');

for t = 2:T;
  x    = find(rand < F(:,S(t-1)));
  S(t) = x(1);
end;

S = S(:); 
