%%
% 
%  PREFORMATTED
%  TEXT
% 
clear all; close all; clc;
% Using function var_Markov_MM.m
% the MM method in Gospodinov and Lkhagvasuren, 2013, for the two
% dimensional VAR(1) process considered in the numerical experiment
% of the  paper.
%function [PN, YN]=var_Markov_MM(A0x,vex,nbar,ntune)
%  A0x      stands for the 2X2 coefficient matrix A.
%  vex      is the 2X2 diagonal matrix, Omega, i.e.  
%               Omega(1,1)=omega_{1,1}^2 and Omega(2,2)=omega_{2,2}^2
%  nbar     is the number of grid points for each i.
%  ntune    is the control variable, where 
%           setting ntune=0 performs the baseline method (MM0), while
%           setting ntune>1 performs the full version of the method,
%           MM. For the examples considered in the paper, ntune was 
%           set to 1000. While higher values of ntune gives a better 
%           approximation, the gain becomes negligible beyond 
%           the value 10000. 

%A0x = [0.8 0.5; 0.1 0.3];
%vex= [0.1 0 ; 0 0.02];
%nbar=9;
%ntune=0;
%[PN, YN]=var_Markov_MM(A0x,vex,nbar,ntune);
%L=4;
%Y=[(1:L)' (11:(10+L))'];
%Y=[(1:L)' (11:(10+L))' (21:20+L)' ];
%Y=[(1:L)' (11:(10+L))' (21:20+L)' (31:(30+L))'];

K=4;  %size(Y,2);
nbar=3;   %size(Y,1);
%A0x_ = random('beta',1,1,K,K); %[0.8 0.5 0 0.001; 0.1 0.3 0 0 ; 0.1 0.2 0.6 0; 0 0.02 0.1  0.02];
%A0x=A0x_*A0x_';
%A0x=diag(random('beta',6,1,K,1));
A0x=[0.9809 0.0028; 0.0410 0.9648];  % el del paper
A0x=[0.9809 0.0028 0.01; 0.0410 0.9648  0.01; 0 -0.01   0.75];  % del paper extendido
A0x=[0.9809 0.0028 0.01 0; 0.0410 0.9648  0.01 0; 0 -0.01  0.75 0.02; 0 0 0.02  0.6];  % del paper aun MAS extendido 
%[eig(A0x_)  eig(A0x) ]
%vex=diag(random('beta',1,1,K,1))/50; %[0.1 0 0 0; 0 0.02 0 0; 0 0 0.0003  0; 0 0 0 0.0001 ];
vex= diag([0.0087,0.0262])^2;  % el del paper
vex= diag([0.0087,0.0262, 0.01])^2;  % el del paper extendido
vex= diag([0.0087,0.0262, 0.01, 0.1])^2;  % el del paper aun MAS extendido
ntune=0;

if K<=2
    [PN_0, YN_0, pmat_0]=var_Markov_MM(A0x,vex,nbar,ntune);
    size(PN_0)
    size(YN_0)
    size(pmat_0)
    fil_max = min(size(PN_0,1),9) ;
    PN_0(1:fil_max,1:fil_max)
    PN_0(end-fil_max+1:end,end-fil_max+1:end)
    %PN_0
end

%[PN, YN]=var_Markov_MM_New(Y,A0x,vex,nbar,ntune);  No es necesario
%guardar pmat, la salvo solo por comparación y chequeos
[PN, YN, pmat]=var_Markov_MM_New(A0x,vex,nbar,ntune);
size(PN)
size(YN)
size(pmat)
fil_max = min(size(PN,1),9) ;
PN(1:fil_max,1:fil_max)
PN(end-fil_max+1:end,end-fil_max+1:end)
%PN

if K<=2
    PN_0-PN
    max(max(PN_0-PN))
    min(min(PN_0-PN))
    [ numel(PN_0)   numel(PN) ]
    [numel(pmat_0)  numel(pmat)]
end

% ntune=1;
% [PN, YN, pmat]=var_Markov_MM_New(A0x,vex,nbar,ntune);
% size(PN)
% size(YN)
% size(pmat)
% fil_max = min(size(PN,1),9) ;
% PN(1:fil_max,1:fil_max)
% PN(end-fil_max+1:end,end-fil_max+1:end)
% %PN
