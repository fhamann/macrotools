function i = getstateindex(s,S) 

%GETSTATEINDEX.M  extrae la posicion del elemento s en el conjunto S 
% 
%  i = getstateindex(s,S) retorna un escalar i, que indica la posicion
%  del vector s, de dimension 1xk, dentro del conjunto de estado S, de
%  dimension nxk, donde k es el numero de variables de estado y n el 
%  numero de valores discretos de las k variables de estado. 
%

%  Code by    : Franz A Hamann ,22-8-13

n = size(S);
s  = ones(n(1),1)*s;
[~,i] = min(abs(s-S));