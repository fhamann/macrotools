function Prd = otpmos2(xr,xd,prob,nb,ns,nz,D,lambda,bzero)

% OTPOS2 Optimal Transition Probability Matrix
%
% Usage:
%           Prd = otpmos2(xr,xd,prob,nb,ns,nz,D,lambda,bzero)
%
%   INPUTS
%       xr     : optimal policy function under repayment (column vector)
%       xd     : optimal policy function under default (column vector)
%       prob   : transition probability matrix for exogenous process
%       nb     : Number of grid points in the first endogenous variable 
%       ns     : Number of grid points in the first endogenous variable 
%       nz     : Number of grid points in the space of exogenous variables
%       D      : Default vector (1 if sovereign defaults)
%       lambda : Exogenous probability of redemption
%       bzero  : Grid position where debt is equal to zero
%
%   OUTPUTS
%       Prd    : Optimal transition probability matrix
%
if length(xr)~=length(xd)
    error('Policy functions must have the same length');
end
xdd=ceil(xd/bzero)*bzero;
n = length(xr);
Pr=zeros(n);
Pd=zeros(n);
PP=[prob;zeros((nb*ns-1)*nz,nz)];
PP2=reshape(PP,nz,n);
PP3=kron(PP2,ones(nb*ns,1));
for i=1:n; 
    Pr(i,xr(i):n)=PP3(i,1:n-xr(i)+1);
    Pd(i,xdd(i):n)=PP3(i,1:n-xdd(i)+1);
end
DD=repmat(D,1,n);
P11=Pr.*(1-DD);
P12=Pd.*DD;
ind=zeros(nb,n);
ind(bzero,:)=1;
Ind=repmat(ind,ns*nz,1);
P21=P11*lambda.*Ind;
P22=Pd*(1-lambda).*Ind;
Prd=[P11 P12; P21 P22];
end