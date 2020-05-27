function [Pstar] = otpm(x,prob,x_order)
% New OTPM Optimal Transition Probability Matrix version 2020
%
% Usage:
%           P = otpm(x,prob,x_order)
%
%   INPUTS
%       x      : optimal policy function nx1 vector
%       prob   : transition probability matrix for exogenous process nexne
%                matrix
%     optional : 
%       x_order: Indicative variable describing how x is sorted. (def =  1)
%              [1]: If x varies first on the exogenous variables and then
%                on the endogenous  
%                2: If x varies first on the endogenous variables and then
%                on the exogenous
%                (Note: Be careful how the gridmake is defined) 
%   OUTPUT
%       Pstar  : Optimal transtion matrix (nxn)

n     = length(x)   ;
nexo  = size(prob,1);
nendo = n/nexo      ;

    if nargin<3
        xt_order = 1;
    else
        xt_order = x_order;
    end

        P1 = repmat(prob,nendo,1);    
    if xt_order ~= 1 
        x = vec(reshape(x,nendo,nexo)');
    end
        Pstar = zeros(n);
    for ix = 1:n
        Pstar(ix,nexo*(x(ix)-1)+1:nexo*x(ix)) = P1(ix,:);
    end
    
    if xt_order ~= 1     
        Pstar = permute(reshape(Pstar,nexo,nendo,nexo,nendo),[2 1 4 3]);
        Pstar = reshape(Pstar,n,n);
    end
