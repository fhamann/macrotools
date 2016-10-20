function P = otpmos(x,prob,lambda,azero,vd,vr,grid,exo,endo)

% OTPMOS Optimal Transition Probability Matrix for Optimal Stopping problem
%
% Usage:
%           P = otpmos(x,pi,lambda,azero,vd,vr,grid,exo,endo)
%
%   INPUTS
%       x      : optimal policy function (column vector)
%       prob   : transition probability matrix for exogenous process
%       lambda : Exogenous probability to re-entry from the outside
%       azero  : grid position for zero debt
%       grid   : grid with the positions of exogenous and endogenous states
%       exo    : column position for the exogenous variable over the grid 
%       end    : column position for the endogenous variable over the grid
%   OUTPUTS
%       P       : Optimal transtion matrix
%
% Note:
%       policy is ordered from exogenous to endogenous variables 

n   = length(x);

P00 = sparse(n,n); P01 = sparse(n,n);
P10 = sparse(n,n); P11 = sparse(n,n);

% For exogenous decisions
for j =1:n
    if (grid(j,endo) == azero)
        for i = 1:n
            P10(i,j) = lambda*prob(grid(i,exo),grid(j,exo));     
            P11(i,j) = (1-lambda)*prob(grid(i,exo),grid(j,exo)); 
        end
    end
end
% For endogenous decisions
for i = 1:n
    for j=1:n
            k=find(grid(:,endo) == x(i));
                if (x(i) == azero && vd(i)>vr(i))
                    P01(i,k) = prob(grid(i,exo),grid(k,exo));
                else
                    P00(i,k) = prob(grid(i,exo),grid(k,exo));
                end
        end
end
 P = [P00,P01;P10,P11];
end