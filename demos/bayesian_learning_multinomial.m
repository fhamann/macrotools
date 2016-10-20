% BAYESIAN_LEARNING  Example of Beta-multinomial model for a Markov chain

% Suppose the exogenous state st takes n discrete values and that the 
% transition probabilities are governed by a Markov transition matrix 
% whose transitions are governed by two probabilities p and q. 
% The states are assumed to be observable, but the probabilities p and q
% are unknown. Our agent learns about them using Bayes's theorem.
%
% Because the states are discrete and take on n values, a
% beta-multiinomial probability model is convenient. We assume the agent 
% has independent beta priors over (p,q) so that f(p,q)=f(p)f(q).
%
% The variable nij represents a counter that records the number of 
% transitions from state t i to j through date t. The parameters nij 
% represent prior beliefs about the frequency 0 of transitions, which may 
% for example come from a training sample. It can be shown that the 
% likelihood function for a batch of data st is proportional to the 
% product of multinomial densities.
%
% Given independent beta priors over p and q and a likelihood function
% that is a product of multinomials, it follows that the posteriors are  
% also independent and have the beta form.

clear all

%% Define the true but unobservable Markov chain
p  = 0.98;       % q_ll                    
q  = 0.98;       % q_hh
m  = 1;         % mean (normalize to 1)
e  = 0.05;      % standard deviation (30% if the mean is 1)
ns = 9;         % number of nodes, usually odd

[X,Q] = markovchain(ns,p,q,e,m)           % Markov chain

%% Simulate the Markov chain

T = 5000;                                   % Sample size
s = simulmarkov(Q,T);                       % Simulate a Markov chain

%% Solve the Bayesian learning problem

% Set the initial counters

n0 = 0.0014;           % n0 -> 0 means that there is no previous knowledge
n = n0*ones(1,ns);     % no previous knowledge for all states
N = repmat(n,ns,1);  

Qt = betamultinomial(s,N);  % Get t=1,...,T transition matrices

%% Compute some moments 

for t=1:T
  EX(t)=X'*ergdist(Qt(:,:,t));     % Anticipated unconditional expectation
  Ec1(t)= condmom(Qt(:,:,t),X,1);  % Conditional expectetaion of state 1
  Ec2(t)= condmom(Qt(:,:,t),X,2);  % Conditional expectetaion of state 2
  Ec3(t)= condmom(Qt(:,:,t),X,3);  % Conditional expectetaion of state 3
  Ec4(t)= condmom(Qt(:,:,t),X,4);  % Conditional expectetaion of state 4
  Ec5(t)= condmom(Qt(:,:,t),X,5);  % Conditional expectetaion of state 5
  Ec6(t)= condmom(Qt(:,:,t),X,6);  % Conditional expectetaion of state 6
  Ec7(t)= condmom(Qt(:,:,t),X,7);  % Conditional expectetaion of state 7
  Ec8(t)= condmom(Qt(:,:,t),X,8);  % Conditional expectetaion of state 8
  Ec9(t)= condmom(Qt(:,:,t),X,9);  % Conditional expectetaion of state 9
end

plot([EX' Ec1' Ec2' Ec3' Ec4' Ec5' Ec6' Ec7' Ec8' Ec9']);