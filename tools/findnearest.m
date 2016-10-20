function [ind,nearest] = findnearest(x,tabvals)

% FINDNEAREST Closest value in a column vector
%
%  Code by    : Franz A Hamann , 3-5-02

for i=1:size(x,2)
    [nearest(i), ind(i)] = min(abs(x(i)-tabvals));
end