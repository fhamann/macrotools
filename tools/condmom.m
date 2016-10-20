function [M,D] = condmom(P,s,i)

% CONDMOM mean and variance of S conditional on a vector of states (s) 
%         with transition probability matrix (P) 
%
% Usage:     [M,D] = condmom(P,s,i)
%
% Inputs:
%
%       P   : n by n transition probability matrix
%       s   : n by 1 column vector of values of the state
%       [i] : optional, i-th position of the state
%
% Output:
%
%       M   : n by 1 vector of conditional means 
%       D   : n by 1 vector of conditional std deviations

M = P*s;
D = sqrt(P*(s.^2)-M.^2);

if nargin>2;
    M = M(i);
    D = D(i);
end 