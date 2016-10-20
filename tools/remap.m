function [x,S] = remap(Sh,S,x)

% REMAP  remaps value/policy function from state-space, S to Sh 
%
% Syntax:          x = remap(Sh,S,x)
%
% Inputs:    Sh : new state-space grid to map to
%            S  : old state-space grid to map from
%            x  : old value/policy function (positions)
%
% Output:    x  : new value/policy function (positions)

i = getindex(Sh,S);
x = x(i);
S = Sh;