clear all;
clc;

T_PX  = 29;
ER1   = 4;
ER2   = 23;

n(1,1) = 0.999999;
n(1,2) = 0.999999;
n(2,1) = 0.000001;  
n(2,2) = 0.000001;

P_X = [0.7 1.3];


px(1:ER1)       = 2;
px(ER1+1:ER2)   = 1;
px(ER2+1:T_PX)  = 2;  


for t=2:T_PX
    
    for i=1:2
       for j=1:2
           if (px(t) == j) && (px(t-1)==i)
               n(i,j,t) = n(i,j,t-1)+1;
           else
               n(i,j,t) = n(i,j,t-1);
           end
       end
    end


end

nhh = squeeze(n(1,1,:));
nhl = squeeze(n(1,2,:)); 
nlh = squeeze(n(2,1,:));
nll = squeeze(n(2,2,:));

E_Fhh = nhh./(nhh+nhl);
E_Fhl = nhl./(nhh+nhl);
E_Fll = nll./(nll+nlh);
E_Flh = nlh./(nll+nlh);


%%%%%%%%%%%%%%%%%%%%%

for t=1:T_PX
EP(t)=P_X*ergdist([ E_Fll(t) 1-E_Fll(t) ; 1-E_Fhh(t) E_Fhh(t)]);
end



VE_Fhh = E_Fhh(1:T_PX);
VE_Fll = E_Fll(1:T_PX);

Fs = [(1:T_PX)' E_Fll  E_Fhh];



save prob_sim_path VE_Fhh VE_Fll T_PX
