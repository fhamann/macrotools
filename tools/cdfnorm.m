function p = cumnorm(x,mu,var)

% CUMNORM Evaluates cumulative probability of univariate normal distribution
%
% Usage
%   p = cumnorm(x,mu,var)
% Input
%   x       : evaluation point
%   mu      : mean of distribution (optional, default=0)
%   var     : variance of distribution (optional, default=1)
% Output
%   p       : cumulative probability

if nargin<2, mu=0; end
if nargin<3, var=1; end

x = (x-mu)/sqrt(var);
z = abs(x)/sqrt(2);
t = 1/(1+0.5*z);

p = 1-0.5*t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t* ...
     (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t* ...
     (1.48851587+t*(-.82215223+t*.17087277)))))))));

if x<0, p=1-p; end
