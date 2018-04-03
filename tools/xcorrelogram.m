function xcf = xcorrelogram(y,x,p)

% XCORRELOGRAM - Cross correlations through p Lags
%
% Usage:
%             xcf = xcorrelogram(y,x,p)
%
% Inputs:
%     y - series to compute acf for, nx1 column vector
%     x - series to compute acf for, nx1 column vector
%     p - total number of lags, scalar
%
% Output:
%     xcf - px1 vector containing  cross correlation function
%

[n1, n2] = size(y) ;
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end

[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
    error('Input p must be scalar and less than length(y)')
end

xcf = zeros(p,1);

N = max(size(y));
ybar = nanmean(y); 
xbar = nanmean(x); 

% Collect ACFs at each lag i
for i = 1:p
   la(i) = lag(y,x,i) ; 
end

yvar = sqrt(nansum((y-ybar).^2)); %yvar = sqrt((y-ybar)'*(y-ybar));
xvar = sqrt(nansum((x-xbar).^2)); %xvar = sqrt((x-xbar)'*(x-xbar));

current = nansum((y-ybar).*(x-xbar))/(yvar*xvar);

for i = 1:p
   le(i) = lead(y,x,i) ; 
end
xcf = [flip(la) current le]';

% ---------------
% SUB FUNCTION
% ---------------
function ta2 = lag(y,x,k)

N = max(size(y)) ;
ybar = nanmean(y); 
xbar = nanmean(x); 
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(x(i-k)-xbar) ;
end

% Denominator, unscaled variance
yvar = sqrt(nansum((y-ybar).^2)); %yvar = sqrt((y-ybar)'*(y-ybar));
xvar = sqrt(nansum((x-xbar).^2)); %xvar = sqrt((x-xbar)'*(x-xbar));

ta2 = nansum(cross_sum) / (yvar*xvar) ;


% ---------------
% SUB FUNCTION
% ---------------
function ta2 = lead(y,x,k)

N = max(size(y)) ;
ybar = nanmean(y); 
xbar = nanmean(x); 
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (x(i)-xbar)*(y(i-k)-ybar) ;
end

% Denominator, unscaled variance
yvar = sqrt(nansum((y-ybar).^2)); %yvar = sqrt((y-ybar)'*(y-ybar));
xvar = sqrt(nansum((x-xbar).^2)); %xvar = sqrt((x-xbar)'*(x-xbar));

ta2 = nansum(cross_sum) / (yvar*xvar) ;