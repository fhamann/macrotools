function [PN, YN, pmat ]=var_Markov_MM(A0x,vex,nbar,ntune)

% This function constructs a finite state Markov chain approximation 
% using the MM method in Gospodinov and Lkhagvasuren, 2013, for the two
% dimensional VAR(1) process considered in the numerical experiment
% of the  paper.
%
% The VAR(1) process:      y' = Ay + epsilon
%
% where var(epsilon) is given by a diagonal matrix Omega.
%
%   INPUT: 
%          A0x stands for the 2X2 coefficient matrix A.
%          vex is the 2X2 diagonal matrix, Omega, i.e.  
%               Omega(1,1)=omega_{1,1}^2 and Omega(2,2)=omega_{2,2}^2
%          nbar  is the number of grid points for each i.
%          ntune is the control variable, where 
%                 setting ntune=0 performs the baseline method (MM0), while
%                 setting ntune>1 performs the full version of the method,
%                 MM. For the examples considered in the paper, ntune was 
%                 set to 1000. While higher values of ntune gives a better 
%                 approximation, the gain becomes negligible beyond 
%                 the value 10000. 
%
%   OUPUT: 
%           PN is the N*-by-N* transition matrix, where N* = nbar^2. The
%                [row k, column j] element is the probability the system 
%                switches from state j to state k. So, the elements of each 
%                column add up to 1.
%           YN is the N*-by-2 matrix of the discrete values of y1 and y2 
%                for N* states. 
%          

if ntune<0
    error('ntune has to be a positive integer (including zero)');
end

if mod(ntune,1)~=0 
    error('ntune has to be a positive integer');;  
end


nx=ntune+1;
n=nbar;
n1=n;
n2=n;



[probtemp, z] = rouwen(0,0,1,n);  
y1=z;
y2=y1;

A0=A0x;

% normalize the initial var so unconditional variances are 1.
[A0new, vynew, vyold, venew]=var_norm(A0, vex);  

vy=vyold;
A=A0new;
ve=venew;


pmat=zeros(2,n,n,n);
px=zeros(2,n,n);

for i=1:n
    for j=1:n
        for k=1:2
            
            mu=A(k,1)*y1(i)+ A(k,2)*y2(j);
            
            vact=ve(k,k);
            
            r=sqrt(1-vact);
            
            [prob1, z] = rouwen(r,0,1,n);
            
            [v1, p, na,nb, dummy_exceed]=cal_mu_fast(mu,vact,n,z);
            
            if nx<2
                
                if na==nb  % if mu is outside of the grids
                    pmat(k,i,j,:)=prob1(:,na);
                else % more relevant case
                    pmat(k,i,j,:)=p*prob1(:,na)+(1-p)*prob1(:,nb);
                end
                                
            else
                
                if na==nb  % if mu is outside of the grids
                    pmat(k,i,j,:)=prob1(:,na);
                else  % begining of the more relevane
                    
                    B=999*ones(nx,6);
                    ixx=0;
                    for ix=1:nx
               vactx=max(0.00000000000001, vact*(1.0-(ix-1)/(nx-1)));
              [v1x, px, nax,nbx, dummy_exceedx]=cal_mu_fast(mu,vactx,n,z);
                        if abs(dummy_exceedx)<0.5
                            ixx=ixx+1;
                            B(ixx,:)=[v1x px nax nbx dummy_exceedx vactx];
                        end
                    end
                    
                    
                    if ixx<1
                        pmat(k,i,j,:)=p*prob1(:,na)+(1-p)*prob1(:,nb);
                    else
                        bvectemp=B(:,1)-vact;
                        dif1=abs(bvectemp);
                        [difx, iz]=min(dif1);
                        
                        
                        pz=B(iz,2);
                        naz=B(iz,3);
                        nbz=B(iz,4);
                        vactz=B(iz,6);
                        
                        rz=sqrt(1-vactz);
                        [probz, z] = rouwen(rz,0,1,n);
                        pmat(k,i,j,:)=pz*probz(:,naz)+(1-pz)*probz(:,nbz);
                        
                    end  % end of the more relevane
                end
            end
            
            
        end
    end
end


% convert the transition probabilities into a conventional form
PN = bigPPP(pmat,n);

ynum=n*n;
ix=0;
Y=zeros(n*n,2);
for i=1:n
    for j=1:n
        ix=ix+1;
        Y(ix,:)=[y1(i) y2(j)];
    end
end

YN=[Y(:,1)*sqrt(vy(1,1)) Y(:,2)*sqrt(vy(2,2))];


