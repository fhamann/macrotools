function pi = ptmn(rho,rhoen,n);

% Creates n point Markov Chain shocks and PTM
% Usage:
%          pi = ptm(rho,rhoen,n)
%
% rho   : autocorrelation of shocks
% rhoen : correlation between shocks
% n     : number of points

aux2 = ((1+rhoen)/(n^2));

if mod(n,2) == 0
    aux1 = flip(linspace(1/n,1,(n^2/2)));
    aux3 = aux1-aux2;
    pi    = [aux2 aux3(:,[2:end]) flip(aux3(:,[2:end])) aux2];
else 
     aux1 = flip(linspace(1/n,1,(n^2+1)/2));
     aux3 = aux1-aux2;
     pi    = [aux2 aux3(:,[2:end]) flip(aux3(:,[2:end-1])) aux2];
end

rho   = rho*ones(n^2,1);
pi    = kron(1-rho,pi);
delta = diag(rho,0);
pi    = pi+delta;

