% SIMULTI   Simulation of the time iteration results
%
% Usage:
%        [hs,hx,hs2x] = simulti(enodes,snodes,sp,x,Q,T,ie,is)
%
% INPUTS
%  enodes  : ne by one nodes for the exogenous state
%  snodes  : ns by one nodes for the endogenous state
%  sp      : ns by ne optimal next period state (obtained from solveti)
%  x       : ns by ne optimal controls (obtained from solveti)
%  Q       : ne by ne Transition Matrix
%  T       : periods of simulation
%  ie      : Initial point on enodes-grid
%  is      : Initial point on snodes-grid
%
% OUTPUTS
%  hs      : Histogram of state variable
%  hx      : Histogram of control variable
%  hs2x    : Histogram of state/control variable
%
% where ne is the number of discretized exogenous states, e. 

<<<<<<< Updated upstream
function [x_t,s_t] = simulti(enodes,snodes,sp,x,Q,T,ie,is)
=======
function [e_t,s_t,sp_t,x_t,ie_t,is_t] = simulti(enodes,snodes,sp,x,Q,T,ie,is)
>>>>>>> Stashed changes

if nargin<6; T = 10000; end

if nargin<7 
    ie = max(1,round(rand*length(enodes)));   % Initial point on e-grid
    is = max(1,round(rand*length(snodes)));   % Initial point on s-grid
end        

e_t  = zeros(T,1);      % exogenous state variable
s_t  = zeros(T,1);      % endogenous state variable
x_t  = zeros(T,1);      % action 
sp_t = zeros(T,1);      % tomorrow's endogenous state
ie_t = zeros(T,1);      % exogenous state variable position
is_t = zeros(T,1);      % endogenous state variable position

rng('default');         % set seed of random number generator 

e = enodes(ie);         % Initial exogenous state variable
s = snodes(is);         % Initial exogenous state variable

CQ = cumsum(Q,2);       % Cumulative (over rows) transition matrix 

for t=1:T
    e_t(t)  = e;        
    s_t(t)  = s;
    x_t(t)  = interp1(snodes,x(:,ie),s);
    sp_t(t) = interp1(snodes,sp(:,ie),s);    
    ie_t(t) = ie;
    s       = sp_t(t);
    is_t(t) = getindex(s,snodes);
    ie      = sum(CQ(ie,:)<rand())+1; % Find next price using a U(0,1)
    e       = enodes(ie); 
end

figure;
hs = histogram(s_t,'Normalization','Probability');
title('Histogram of state variable')

figure;
hx = histogram(x_t,'Normalization','Probability');
title('Histogram of control variable')

figure;
hs2x = histogram(s_t./x_t,'Normalization','Probability');
title('Histogram of state/control variable')
