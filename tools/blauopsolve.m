function [] = blauopsolve(u,disc,p,s,P,S,path,n0,s0)

[n,m] = size(pie);

T   = length(path);
p_t = p(path)';            % maps position i to vector of prices pt (Tx1)

% Declare vector and matrices for storage for t=1,...,T
H   = zeros(n,n,T);   % n x n x T optimal TPMs for for t=1,...,T
ix  = zeros(n,T);     % each column is vector of policy index for t=1,...,T
x_t = zeros(T,1);     % each column is vector of extraction for t=1,...,T
s_t = zeros(T+1,1);   % each column is vector of reserves for t=1,...,T
prt = zeros(T,1);


s_t(1) = s0;                       % initial stock of reserves

Qt = betabinomial(path,n0);            % Q_t TPM's for t=1,...,T

i = getindex([p_t(1) s_t(1)],[P S]);   % initial state index

for t=1:T
    iold = i;
    Qbig = kron(speye(ns,ns),repmat(Qt(:,:,t),ns,1));
    [~,ix(:,t),H(:,:,t)] = solvedp(u,Qbig,disc,'policy');
    s_t(t+1) = s(ix(iold,t));
    x_t(t) = max(0,s_t(t) + d - s_t(t+1));
    prt(t) = p_t(t).*x_t(t)-kappa*(x_t(t).^2)./s_t(t);
    i      = getindex([p_t(t+1) s_t(t+1)],[P S]);     % update state index
