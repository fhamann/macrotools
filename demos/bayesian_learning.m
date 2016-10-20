% BAYESIAN_LEARNING  Example of Beta-binomial model for a Markov chain

% Suppose the exogenous state st takes two values (0 and 1) and that the 
% transition probabilities are governed by a Markov transition matrix 
% P = [p 1-p; 1-q q]. The states are assumed to be observable, but the 
% transition probabilities are unknown. Our agent learns about them using 
% Bayes?s theorem.
%
% Because the states are discrete and take on two values, a
% beta-binomial probability model is convenient. We assume the agent has 
% independent beta priors over (p,q) so that f(p,q)=f(p)f(q) where
%
%             f(p) ~ p^(n11-1)*(1-p)^(n10-1)
%             f(q) ~ q^(n00-1)*(1-q)^n(01-1)
%
% The variable nij represents a counter that records the number of 
% transitions from state t i to j through date t. The parameters nij 
% represent prior beliefs about the frequency 0 of transitions, which may 
% for example come from a training sample. It can be shown that the 
% likelihood function for a batch of data st is proportional to the 
% product of binomial densities.
%
% Given independent beta priors over p and q and a likelihood function
% that is a product of binomials, it follows that the posteriors are also 
% independent and have the beta form.


%% Define the true, but unobservable Markov chain

p = 0.94;       % prob of moving up after low                    
q = 0.98;       % prob of moving down after low
m = 1;          % mean (normalize to 1)
e = 0.5;        % standard deviation (30% if the mean is 1)

[X,Q] = markovchain(2,p,q,e,m)              % Markov chain

%% Simulate the Markov chain

T = 5000;                                    % Sample size
s = simulmarkov(Q,T);                         % Simulate a Markov chain

%% Solve the Bayesian learning problem

% Set the initial counters

%n  = [0.0014 0.0014; 0.0014 0.0014];  % assume no previous knowledge
n  = Q;                               % assume perfect previous knowledge 
  
Qt = betamultinomial(s,n);              % Get t=1,...,T transition matrices

%% Compute some moments 

for t=1:T
  EX(t)=X'*ergdist(Qt(:,:,t));     % Anticipated unconditional expectation
  EcL(t)= condmom(Qt(:,:,t),X,1);  % Conditional expectetaion of low state
  EcH(t)= condmom(Qt(:,:,t),X,2);  % Conditional expectetaion of low state
end

figure(1); plot([EX' EcL' EcH']);  % Plot the anticipated expectation