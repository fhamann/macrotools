% EJEMPLO2 - BROCK-MIRMAN MODEL FOR TIME ITERATION - Inverse
function [out1,out2,out3] = ejemplo2(flag,k,c,e,mu,kp,alpha);
switch flag
case 'f'  % REWARD FUNCTION
   out1 = log(c);                         % f
   out2 = (c).^(-1);                      % fx
   out3 = zeros(size(k));                 % fs
case 'g'  % STATE TRANSITION FUNCTION
   out1 = e.*(k.^alpha)-c;                % g
   out2 = -ones(size(k));                 % gx
   out3 = alpha*e.*(k.^(alpha-1));        % gs
case 'i'
   out1 = mu.^(-1);                        %fxi
   out2 = e.*(k.^alpha)-kp;  
end 