function varargout = crosscorr(y1,y2,numLags,numSTD)
%CROSSCORR Sample cross-correlation
%
% Syntax:
%
%   [xcf,lags,bounds] = crosscorr(y1,y2)
%   [xcf,lags,bounds] = crosscorr(y1,y2,numLags,numSTD)
%   crosscorr(...)
%
% Description:
%
%   Compute the sample cross-correlation function (XCF) between univariate,
%   stochastic time series. When called with no output arguments, CROSSCORR
%   plots the XCF sequence with confidence bounds.
%
% Input Arguments:
%
%   y1 - Vector of observations of the first univariate time series for 
%     which the sample XCF is computed or plotted. The last element of y1 
%     contains the most recent observation.
%
%   y2 - Vector of observations of the second univariate time series for 
%     which the sample XCF is computed or plotted. The last element of y2 
%     contains the most recent observation.
%
% Optional Input Arguments:
%
%   numLags - Positive integer indicating the number of lags of the XCF to
%     compute. If empty or missing, the default is to compute the XCF at
%     lags 0, +/-1, +/-2,...,+/-T, where T is the minimum of 20 or one less
%     than the length of the shortest series.
%
%   numSTD - Positive scalar indicating the number of standard deviations
%     of the sample XCF estimation error to compute assuming y1/y2 are 
%     uncorrelated. If empty or missing, the default is numSTD = 2
%     (approximate 95% confidence).
%
% Output Arguments:
%
%   xcf - Sample cross-correlation function between y1 and y2. xcf is a
%     vector of length 2*numLags+1 corresponding to lags 0, +/-1, +/-2, ...
%     +/-numLags. The center element of xcf contains the zeroth-lag cross-
%     correlation. xcf will be a row (column) vector if y1 is a row 
%     (column) vector.
%
%   lags - Vector of lags corresponding to xcf (-numLags to +numLags).
%
%   bounds - Two-element vector indicating the approximate upper and lower
%     confidence bounds, assuming the input series are uncorrelated.
%
% Example:
%
%   % Create a random sequence of 100 Gaussian deviates and a delayed
%   % version lagged by 4 samples. Observe the XCF peak at the 4th lag:
%
%     x = randn(100,1);    % 100 Gaussian deviates ~ N(0,1)
%     y = lagmatrix(x,4);  % Delay it by 4 samples
%     y(isnan(y)) = 0;     % Replace NaN's with zeros
%     crosscorr(x,y)       % It should peak at the 4th lag
%
% Reference:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
% See also AUTOCORR, PARCORR, FILTER.

%   Copyright 1999-2010 The MathWorks, Inc.   

if nargin < 2
    
   error(message('econ:crosscorr:UnspecifiedInput'))
     
end

% Ensure the sample data are vectors:

[rows,columns] = size(y1);

if ((rows ~= 1) && (columns ~= 1)) || (rows*columns < 2)
    
   error(message('econ:crosscorr:NonVectorSeries1'))
     
end

[rows,columns] = size(y2);

if ((rows ~= 1) && (columns ~= 1)) || (rows*columns < 2)
    
   error(message('econ:crosscorr:NonVectorSeries2'))
     
else

end

rowSeries = (size(y1,1) == 1);

y1 = y1(:); % Ensure a column vector
y2 = y2(:); % Ensure a column vector

N = min(length(y1),length(y2)); % Sample size
defaultLags = 20; % Recommendation of [1]

% Ensure numLags is a positive integer or set default:

if (nargin >= 3) && ~isempty(numLags)
    
   if numel(numLags) > 1
       
      error(message('econ:crosscorr:NonScalarLags'))
        
   end
   
   if (round(numLags) ~= numLags) || (numLags <= 0)
       
      error(message('econ:crosscorr:NonPositiveInteger'))
        
   end
   
   if numLags > (N-1)
       
      error(message('econ:crosscorr:InputTooLarge'))
        
   end
   
else
    
   numLags = min(defaultLags,N-1); % Default
   
end

% Ensure numSTD is a positive scalar or set default:

if (nargin >= 4) && ~isempty(numSTD)
    
   if numel(numSTD) > 1
       
      error(message('econ:crosscorr:NonScalarSTDs'))
        
   end
   
   if numSTD < 0
       
      error(message('econ:crosscorr:NegativeSTDs'))
        
   end
   
else
    
   numSTD = 2; % Default
   
end

% The FILTER command could be used to compute the XCF, but FFT-based 
% computation is significantly faster for large data sets.

y1 = y1-mean(y1);
y2 = y2-mean(y2);
L1 = length(y1);
L2 = length(y2);

if L1 > L2
    
   y2(L1) = 0;
   
elseif L1 < L2
    
   y1(L2) = 0;
   
end

nFFT = 2^(nextpow2(max([L1 L2]))+1);
F = fft([y1(:) y2(:)],nFFT); 

ACF1 = ifft(F(:,1).*conj(F(:,1)));
ACF2 = ifft(F(:,2).*conj(F(:,2)));

xcf = ifft(F(:,1).*conj(F(:,2)));
xcf = xcf([(numLags+1:-1:1) (nFFT:-1:(nFFT-numLags+1))]);
xcf = real(xcf)/(sqrt(ACF1(1))*sqrt(ACF2(1)));

lags = (-numLags:numLags)';
bounds = [numSTD;-numSTD]/sqrt(N);


if nargout == 0

%  Plot the sample XCF:

   lineHandles = stem(lags,xcf,'filled','r-o');
   set(lineHandles(1),'MarkerSize',4)
   grid('on')
   xlabel('Lag')
   ylabel('Sample Cross Correlation')
   title('Sample Cross Correlation Function')
   hold('on')
   
%  Plot the confidence bounds (horizontal lines) under the hypothesis that
%  the underlying series are uncorrelated.

   a = axis;

   plot([a(1) a(1); a(2) a(2)],[bounds([1 1]) bounds([2 2])],'-b');

   plot([a(1) a(2)],[0 0],'-k');
   hold('off')

else

%  Re-format outputs for compatibility with the y1 input. When y1 is input
%  as a row vector, then pass the outputs as a row vectors; when y1 is a
%  column vector, then pass the outputs as a column vectors.

   if rowSeries
       
      xcf = xcf';
      lags = lags';
      bounds = bounds';
      
   end

   varargout = {xcf,lags,bounds};

end