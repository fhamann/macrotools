function  dist_stats = dstats(x,DSetup)
% dstats function to compute distributional statistics 
% USAGE
%           dist_stats = dstats(x,DSetup)
% Inputs
%   x: n by 2 array. Contains the n sorted values of the variable in the first column and the n values of the sorted ergodic distribution in the second
%   Optional:
%   DSetup: nset by 2 array. Contains nset intervals of the distribution
%           If no information is provided, the predetermined values are: 
%           bottom 40, 40 to 60, 60 to 80,  top 20, top 10, top 5, top 1              
%           The DSetup syntaxis of the predetermined is given by:
%                DSetup =   [0.00 0.40;
%                            0.40 0.60;
%                            0.60 0.80;
%                            0.80 1.00;
%                            0.90 1.00;
%                            0.95 1.00;
%                            0.99 1.00];
% Output 
%  dist_stats: nset by 1 array. Contains the distributional statistics.
% 
% 
% Written by F. Hamann & Juan Méndez (2020). 
% Feel free to copy, change and distribute.
% 

%% Preamble Woodwork
  n = (size(x,1))   ; % Space size
  x = [x,zeros(n,3)]; % Create a matrix to store the cumulative distributions
%% Compute cumulative distributions  
  x(:,3) = cumsum(x(:,2))           ; % Cumulative Ergodic distribution
  x(:,4) = cumsum(x(:,1).*x(:,2))   ; % Cumulative distribution of x
  x(:,5) = x(:,4)/x(n,4)            ; % Normalization to 1 with the Upper value  
%% Compute ditributional statistics
    if nargin<2 % % Use a Predetermined Setup: bottom 40, 40 to 60, 60 to 80,  top 20, top 10, top 5, top 1              
        DSetup = [0.00 0.40;0.40 0.60;0.60 0.80;0.80 1.00;0.90 1.00;0.95 1.00;0.99 1.00];
    end
    nset       = size(DSetup,1)     ;  % Size of the desired distributional statistics
    dist_stats = zeros(nset,1)      ;  % Storing Matrix
    for iset = 1:nset
%       Find in the ergodic distribution the lower and upper value of the distributional statistic        
        [~,a] = min(abs(DSetup(iset,1)-x(:,3))); 
        [~,b] = min(abs(DSetup(iset,2)-x(:,3)));
%       Then evaluate those positions in the normalized cumulative distribution of the variable            
        dist_stats(iset,1) =  x(b,5)-x(a,5);    
    end
    
    
    
