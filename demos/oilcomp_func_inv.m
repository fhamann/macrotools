function [out1,out2,out3] = oilcomp(flag,s,x,e,mr,sp,phi,gama,d);
 switch flag
  case 'f'
     out1 = e.*x - (phi/2)*((x./s).^gama).*x;                    % f
     out2 = e - (phi/2)*((gama+1)*(x./s).^(gama));               % fx
     out3 = (phi/2)*(gama)*((x./s).^(gama+1));                  % fs  
  case 'g'
     out1 = s + d - x;      % g
     out2 = -ones(size(s)); % gx
     out3 = ones(size(s));  % gs 
  case 'i'
     out1 = (((e - mr).*(s.^gama))/((phi/2)*(gama+1))).^(1/gama);
     out2 = s+d-sp;
 end