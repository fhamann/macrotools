function [y,P,d]=tauchen(n,mu,A,S,m)
% N-Dimensional Tauchen method to approximate VAR(1) process by Markov chain
% y(t+1) = (I-A)*mu + A*y(t) + e(t),  where e(t)~Normal(0,S)
% where y, mu, e are M x 1 vectors, and A is M x M matrix
% Input Arguments:
%      n - M x 1 vector of number of states in Markov chain components
%     mu - M x 1 vector of unconditional mean of process
%      A - M x M autocorrelation coefficient matrix
%      S - M x 1 vector of standard deviations of innovations
% Optional Input Argument:
%      m - M x 1 vector of bandwidths (# of std.dev. +- mean) (default=3)
% Output Arguments:
%      y - M x N Markov chain, N=prod(n)
%      P - N x N Markov transition matrix (tomorrow x today)
% Optional Output Argument:
%      d - 1 x N Stationary distribution of Markov chain
% Requires Control System Toolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Iskander Karibzhanov
%          Department of Economics
%          University of Minnesota
%          karib003@umn.edu
K=cumprod(n);
M=numel(n); % # of variables in VAR
N=K(M); % # of states in Markov chain
if nargin<5; m=3*ones(M,1); end
s=m(:).*sqrt(diag(dlyap(A,diag(S.^2)))); % std.dev. of y
z=cell(M,1); y=nan(M,N);
for i=1:M
    if n(i)>1, z{i}=linspace(-1,1,n(i))*s(i); else z{i}=0; end
    y(i,:)=reshape(repmat(z{i},K(i)/n(i),N/K(i)),1,N);
end
s=s./(n(:)-1); P=1;
for i=1:M
    h=normcdf(bsxfun(@minus,z{i}'+s(i),A(i,:)*y)/S(i));
    h(n(i),:)=1; h=permute([h(1,:);diff(h,1,1)],[3 1 2]);
    P=reshape(repmat(h,K(i)/n(i),N/K(i)),N,N).*P;
end
y=bsxfun(@plus,y,mu(:));
if nargout>2
    [d,~]=eigs(P,1,1);
    if sum(d)<0; d=-d; end; d(d<0)=0; d=d'/sum(d);
end