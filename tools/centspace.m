function x = centspace(xm,epsilon,n)

% CENTSPACE  Centered linear simetric space around a value

x = xm+linspace(-epsilon,epsilon,n)';
