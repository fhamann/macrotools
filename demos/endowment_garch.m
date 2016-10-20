%% ENDOWMENT.M  Small Open GARCH-Endowment Economy with Incomplete Markets

% Written by F. Hamann. Feel free to copy, change and distribute.

 clear all

%% Parametros
 sigma	= 2;         % coeficiente de aversion al riesgo 2
 beta	= 0.985;     % factor de descuento
 R      = 1.0125;    % tasa de retorno del activo (vs. R=1)

if R*beta>=1; display('Set beta*R<1 for convergence'); break; end;

%% Cadena de Markov de y
 m  = 3;
 n  = 3;
                                      % Benchmark parameter values:
 b0     = 0.0015;                     % b0     = 0.0015;
 b1     = 0.5;                        % b1     = 0.5;
 b2     = 0.3;                        % b2     = 0.3;
 theta  = 0.3;                        % theta  = 0.3;
 lambda = 0.2;                        % lambda = 0.2;
 h_1    = b0/(1-b1-b2*(1+(theta)^2)); % h_1    = b0/(1-b1-b2*(1+(theta)^2))
 p_0    = 1;                          % p_0    = 1;
 I_p    = 0.35;                        % I_p    = 0.2882;
 q_1    = log(h_1);                   % q_1    = log(h_1);
 I_q    = 5;                          % I_q    = 5.7167;
 
 [PI,y,q] = garch(b0,b1,b2,theta,lambda,p_0,I_p,q_1,I_q,m,n)

 ymean = y'*markov(PI);
 ny    = length(y);

return
%% Grilla de B = {bmin<b2<b3<...<bmax}
 nb    = 500;
 bmin  = -15;
 bmax  =  0;
 b     = linspace(bmin,bmax,nb)';

%% Estado-espacio S = YxB
 [Y,B] = gridmake(y,b); [Q,B] = gridmake(q,b);
 bb    = b;

%% Conjunto de consumo factible C>=0 y funcion de utilidad
 C     = zeros(nb*ny,nb);

 for i=1:nb    
    C(:,i)=Y+R*B-bb(i);   % Y no-almacenable, B bono no-contingente
 end
 u  = (C.^(1-sigma)-1)./(1-sigma); u(C<=0)=-Inf; clear C;

%% Matriz de probabilidad de transicion
 P = kron(speye(nb,nb),repmat(PI,nb,1));

%% Ecuacion de Bellman
 [v,x,pstar] = solvedp(u,P,beta,'value'); 
 bnext       = b(x);
 clear P u;

%% Calculo de la distribucion estacionaria
 pie   = ergdist(pstar);
 F     = cumsum(pie);

 cons  = Y+R*B-bnext;

 bmean = pie'*bnext;
 cmean = pie'*cons;

%% Propiedades del modelo

 %plotdp(v,x,pstar,Y,B,y,b); return


%% Simulacion
 T     = 500;      
 s0    = findnearest(bmean,B);   
 spath = simulmarkov(pstar,T,s0);

 ypath = Y(spath);
 cpath = ypath + R*B(spath)-b(x(spath));
 qpath = Q(spath);

 sd.y = std(ypath);
 sd.c = std(cpath);

 figure(1); plot([ypath cpath])
 figure(2); plot(B,pie)

%% Estadisticas del modelo
 display(' ')
 display('Model Summary ')

 rowname = ['y mean       ';
            'c mean       ';
            'b mean       ';
            'b/y          ';
            'Y stdv       ';
            'C stdv       ';
            'Csdv/Ysdv    '];
         
 values=[ymean;cmean;bmean;bmean/ymean;sd.y;sd.c;sd.c/sd.y];

 [rowname num2str(values)]