function plotdp(v,x,pstar,E,A,e,a)

% PLOTDP Plots solution to DP problem
%
% Syntax: plotdp(v,x,pstar,E,A,e,a)
%
% Inputs:  v     :  value function (from solvedp) 
%          x     :  policy function (from solvedp) 
%          pstar :  optimal transition matrix from solvedp 
%          E,A   :  state-space grid 
%          e,a   :  values for exogenous and endogenous states 
% 
% Outputs: a 4x4 plot of value & policy functions and stationary
%          density & distribution of the state


ne = length(e);
na = length(a);

amin = min(a);
amax = max(a);

ax = a(x);
pi    = ergdist(pstar);
F     = cumsum(pi);

pir    = reshape(pi,ne,na)';
axr    = reshape(ax,ne,na)';
vr     = reshape(v,ne,na)';

figure(1)
subplot(221);h=plot(a,axr);
set(h,'color','k','linewidth',1);
axis([amin amax amin amax])
xlim=get(gca,'xlim'); ylim=get(gca,'ylim');
h=line(xlim,ylim);
set(h,'color','k','linestyle','--');
title('Policy function')

subplot(222);
plot(a,vr);
xlabel('a'); ylabel('v(a,e)')
title('Value function')

subplot(223)
plot(A,pi)
title('Stationary density')

subplot(224)
plot(A,F)
title('Stationary distribution')