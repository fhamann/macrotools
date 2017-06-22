% EJEMPLO - HUGGETT PARTIAL EQUILIBRIUM MODEL FOR TIME ITERATION
function [out1,out2,out3] = ejemplo(flag,b,c,e,gamma,q);
switch flag
case 'f'  % REWARD FUNCTION
   out1 = ((c).^(1-gamma))/(1-gamma);     % f
   out2 = (c).^(-gamma);                  % fx
   out3 = zeros(size(b));                 % fs
case 'g'  % STATE TRANSITION FUNCTION
   out1 = (1/q)*(b+e-c);                  % g
   out2 = -(1/q)*ones(size(b));           % gx
   out3 = (1/q)*ones(size(b));            % gs
end 