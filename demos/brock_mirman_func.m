% EJEMPLO2 - BROCK-MIRMAN MODEL FOR TIME ITERATION
function [out1,out2,out3] = ejemplo2(flag,s,x,e,alpha);
switch flag
case 'f'  % REWARD FUNCTION
   out1 = log(x);                         % f
   out2 = (x).^(-1);                      % fx
   out3 = zeros(size(s));                 % fs
case 'g'  % STATE TRANSITION FUNCTION
   out1 = e.*(s.^alpha)-x;                % g
   out2 = -ones(size(s));                 % gx
   out3 = alpha*e.*(s.^(alpha-1));        % gs
end 