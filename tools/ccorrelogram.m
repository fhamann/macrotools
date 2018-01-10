function ta = ccorrelogram(y,x,p)
% ACF - Compute Autocorrelations Through p Lags
% >> myacf = acf(y,p) 
%
% Inputs:
% y - series to compute acf for, nx1 column vector
% p - total number of lags, 1x1 integer
%
% Output:
% myacf - px1 vector containing autocorrelations
%        (First lag computed is lag 1. Lag 0 not computed)
%
%
% A bar graph of the autocorrelations is also produced, with
% rejection region bands for testing individual autocorrelations = 0.
%
% Note that lag 0 autocorelation is not computed, 
% and is not shown on this graph.
%
% Example:
% >> acf(randn(100,1), 10)
%


% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(y) ;
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end

[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end



% -------------
% BEGIN CODE
% -------------

ta = zeros(p,1) ;
global N 
N = max(size(y)) ;
global ybar 
ybar = mean(y); 
global xbar 
xbar = mean(x); 

% Collect ACFs at each lag i
for i = 1:p
   la(i) = lag(y,x,i) ; 
end

yvar = sqrt((y-ybar)'*(y-ybar) );
xvar = sqrt((x-xbar)'*(x-xbar) );
current = sum((y-ybar).*(x-xbar))/(yvar*xvar);

for i = 1:p
   le(i) = lead(y,x,i) ; 
end
ta = [flip(la) current le]';



% ---------------
% SUB FUNCTION
% ---------------
function ta2 = lag(y,x,k)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
global ybar
global xbar
global N
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(x(i-k)-xbar) ;
end

% Denominator, unscaled variance
yvar = sqrt((y-ybar)'*(y-ybar) );
xvar = sqrt((x-xbar)'*(x-xbar) );

ta2 = sum(cross_sum) / (yvar*xvar) ;

% ---------------
% SUB FUNCTION
% ---------------
function ta2 = lead(y,x,k)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
global ybar
global xbar
global N
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (x(i)-xbar)*(y(i-k)-ybar) ;
end

% Denominator, unscaled variance
yvar = sqrt((y-ybar)'*(y-ybar) );
xvar = sqrt((x-xbar)'*(x-xbar) );

ta2 = sum(cross_sum) / (yvar*xvar) ;


