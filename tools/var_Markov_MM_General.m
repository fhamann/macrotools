% This function constructs a finite state Markov chain approximation 
% using the MM method in Gospodinov and Lkhagvasuren, 2013, for the two
% dimensional VAR(1) process considered in the numerical experiment
% of the  paper.
%  https://sites.google.com/site/dlkhagva/var_mmm
%  http://qed.econ.queensu.ca/jae/datasets/gospodinov001/
% The VAR(1) process: y'=Ay+epsilon,
% where var(epsilon) is given by a daigonal matrix Omega.
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
%           type  =1 for all version =2 for new version with posible
%                 diferent number of states by each variable
%                 if type =2 TransititonProbabiliyMatrix_General_2.m is used
%              type=0 use TransititonProbabiliyGeneralMatrix.m  older version
%   OUPUT: 
%           PN is the N*-by-N* transition matrix, where N* = nbar^2. The
%                [row k, column j] element is the probability the system 
%                switches from state j to state k. So, the elements of each 
%                column add up to 1.
%           YN  is the N*-by-2 matrix of the discrete values of y1 and y2 
%                for N* states. 
%          Pmat  a tensor with intermedaite results previous to Pmat
%
function [PN, YN, Pmat]=var_Markov_MM_General(A0x,vex,nbar,ntune, type)
K=size(A0x,1);
if ntune<0
    error('ntune has to be a positive integer (including zero)');
end
if mod(ntune,1)~=0 
    error('ntune has to be a positive integer');
end
if length(nbar)==1
    %nx=ntune+1; n=nbar;
    n=ones(K,1)*nbar;
else
    if length(nbar) ~= K
       error('length(nbar) ~= dim(A)');
    end
    n=nbar;
end
%n1=n;
%n2=n;
y_i={};
%n
for i=1:K
  %[~, z] = rouwen(0,0,1,n); 
  %n(i)
  [~, z] = rouwen(0,0,1,n(i));
  %y1=z; y2=y1;
  %y_i = kron(z,ones(1,K)) ;
  y_i{i} =  z ;
end
  %y_i
  A0=A0x;

% normalize the initial var so unconditional variances are 1.
[A0new, ~, vyold, venew]=var_norm(A0, vex);  
%A0new
vy=vyold;
A=A0new;
ve=venew;
%ynum=n*n;
%ix=0;
%y_i
%Y=Permutacion_columnas(y_i);
Y=Permutacion_columnas_General(y_i);
%Y
%function P_out = TransititonProbabiliyGeneralMatrix(A,ve,n,ntune)
%Pmat = TransititonProbabiliyGeneralMatrix(A,ve,n,ntune,y_i,Y);
if type == 2
    Pmat = TransititonProbabiliyMatrix_General_2(A,ve,n,ntune,y_i,Y);
elseif type == 1
    Pmat = TransititonProbabiliyMatrix_General(A,ve,n,ntune,y_i,Y);
elseif type == 0
    Pmat = TransititonProbabiliyGeneralMatrix(A,ve,n(1),ntune,y_i,Y);
else
    error('You must choice type=0, type=1  or type=2')
end
%size(Pmat)
% convert the transition probabilities into a conventional form
%PN = bigPPP(pmat,n);
%PN = bigPPP_New(Pmat,n);
PN = bigPPP_General(Pmat,n);
% estos condicionales son necesarios cuando el # de estados para cada
% variable no es una constante
if max(all(PN==0,1)) ~= 0
    'AJA COLUMNAS'
    PN(:,all(PN==0,1))=[];
end
if max(all(PN==0,2)) ~= 0
    'AJA FILAS'
    PN(all(PN==0,2),:)=[];
end
%PN

YN=Y(:,1)*sqrt(vy(1,1));
for i=2:K
  YN=[YN Y(:,i)*sqrt(vy(i,i))];
end

PN = PN';