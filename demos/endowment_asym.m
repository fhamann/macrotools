% ENDOWMENT_ASYM.M  Small Open Endowment Economy with Incomplete Markets
%                   with large downside risk

% Written by F. Hamann. Feel free to copy, change and distribute.

 clear all

%% Parametros
 sigma	= 2;         % coeficiente de aversion al riesgo 2
 beta	= 0.985;     % factor de descuento
 R      = 1.0125;    % tasa de retorno del activo (vs. R=1)

 if R*beta>=1; display('Set beta*R<1 for convergence'); break; end;

%% Cadena de Markov de y
 ny = 3;
 p1 = 0.8;
 p2 = 0.85;
 q  = 0.8;
 e  = 0.2;
 ym = 1;
 a  = 1;
 
 [y1,PI1] = markovchain(ny,p1,q,e,ym)
 [y2,PI2] = markovchain(ny,p2,q,e,ym,a)
  
%% Grilla de B = {bmin<b2<b3<...<bmax}
 nb    = 500;
 bmin  = -15;
 bmax  =  5;
 b     = linspace(bmin,bmax,nb)';

%% Estado-espacio S = YxB
 [Y1,B] = gridmake(y1,b);
 [Y2,B] = gridmake(y2,b);

%% Conjunto de consumo factible C>=0 y funcion de utilidad
 C     = zeros(nb*ny,nb);

 for i=1:nb    
    C(:,i)=Y1+R*B-b(i);   % Y no-almacenable, B bono no-contingente
 end
 u1  = (C.^(1-sigma)-1)./(1-sigma); u1(C<=0)=-Inf; 
 
 for i=1:nb    
    C(:,i)=Y2+R*B-b(i);   % Y no-almacenable, B bono no-contingente
 end
 u2  = (C.^(1-sigma)-1)./(1-sigma); u2(C<=0)=-Inf; 
  
 clear C;

%% Matriz de probabilidad de transicion
 P1 = kron(speye(nb,nb),repmat(PI1,nb,1));
 P2 = kron(speye(nb,nb),repmat(PI2,nb,1));

%% Ecuacion de Bellman
 [v1,x1,pstar1] = solvedp(u1,P1,beta,'policy');  
 [v2,x2,pstar2] = solvedp(u2,P2,beta,'policy'); 
 
 clear P1 P2 u1 u2;

%% Calculo de la distribucion estacionaria
 ymean1 = y1'*markov(PI1);
 ysdv1  = sqrt(markov(PI1)'*(y1-ymean1).^2);    
 
 ymean2 = y2'*markov(PI2);
 ysdv2  = sqrt(markov(PI2)'*(y2-ymean2).^2);                       

 pie1   = ergdist(pstar1);
 F1     = cumsum(pie1);

 cons1  = Y1+R*B-b(x1);
 CA1    = b(x1)-B;

 bmean1 = pie1'*b(x1);
 bsdv1  = sqrt(pie1'*(b(x1)-bmean1).^2);                       

 cmean1 = pie1'*cons1;
 csdv1  = sqrt(pie1'*(cons1-cmean1).^2);                       

 CAmean1  = pie1'*CA1;
 CAsdv1   = sqrt(pie1'*(CA1-CAmean1).^2);

 pie2   = ergdist(pstar2);
 F2     = cumsum(pie2);

 cons2  = Y2+R*B-b(x2);
 CA2    = b(x2)-B;

 bmean2 = pie2'*b(x2);
 bsdv2  = sqrt(pie2'*(b(x2)-bmean2).^2);                       

 cmean2 = pie2'*cons2; 
 csdv2  = sqrt(pie2'*(cons2-cmean2).^2);                       

 CAmean2  = pie2'*CA2;
 CAsdv2   = sqrt(pie2'*(CA2-CAmean2).^2);

% Plots stationary distributions  
 figure(1); plot(B,[F1 F2]);     legend('symmetric','asymmetric')
 figure(2); plot(B,[pie1 pie2]); legend('symmetric','asymmetric')


%% Simulation

s1 = simulmarkov(PI1,500,2);
s2 = simulmarkov(PI2,500,2);
figure(3); plot([Y1(s1) Y2(s2) cons1(s1) cons2(s2)]);

%% Table report
names   = ['Output               ';
           'Consumption          ';    
           'Current account      ';
           'Net foreign assets   ';
           'NFA to GDP           '];

means1 = [ymean1;cmean1;CAmean1;bmean1;bmean1/ymean1];
means2 = [ymean2;cmean2;CAmean2;bmean2;bmean2/ymean2];

espacio = ['   ';'   ';'   ';'   ';'   '];
stddevs1 = [ysdv1;csdv1;CAsdv1;bsdv1; NaN];
stddevs2 = [ysdv2;csdv2;CAsdv2;bsdv2; NaN];

display('Means & Std.Devs of Standard vs Stoch. Volatility Model')
display('Means in cols 1 & 2; Std.Devs in cols 3 & 4')

TABLE = [names num2str(means1)   espacio num2str(means2)  espacio...
               num2str(stddevs1) espacio num2str(stddevs2) ]
