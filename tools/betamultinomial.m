function Ps  = betamultinomial(s,n)

% BETABINOMIAL  Beta-binomial model for a n-state Markov chain 
%
% Usage:      Ps  = betamultinomial(s,n)
%
% Input:
%
%       s  : T by 1 column vector of values of the states
%       n  : n by n matrix of initial counters 
%
% Output:
%
%       Ps : sequence t=1,...,T of transition probability matrices as a
%            function of the vector s conditional on past history up to t


T  = length(s);
ns = length(n);

Ps = zeros(ns,ns,T);
Ps1 = zeros(ns,ns,T);

for t=2:T   
    for i=1:ns
       for j=1:ns
           if (s(t) == j) && (s(t-1)==i)
               n(i,j,t) = n(i,j,t-1)+1;
           else
               n(i,j,t) = n(i,j,t-1);
           end
       end
    end
end

for i=1:ns
   for j=1:ns
      eval(['n_',int2str(i),int2str(j),'=squeeze(n(',int2str(i),',',int2str(j),',:));']);
      if j==1
          eval(['n_',int2str(i),'_tot=n_',int2str(i),int2str(j),';']);
      elseif j>1
          eval(['n_',int2str(i),'_tot=n_',int2str(i),'_tot+n_',int2str(i),int2str(j),';']);
      end
   end
   for k=1:ns
      eval(['Ep_',int2str(i),int2str(k),'=n_',int2str(i),int2str(k),'./(n_',int2str(i),'_tot);']);
   end
end

for i=1:ns
   for j=1:ns
      eval(['Ps(',int2str(i),',',int2str(j),',:)=Ep_',int2str(i),int2str(j),';']);
   end
end
end