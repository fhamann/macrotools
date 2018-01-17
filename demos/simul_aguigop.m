% AGUIAR-GOPINATH SIMULATION
%
% Simulates stochastic oil price and endowment & sovereign default model
% Sovereign owns and operates oil extraction company decides optimal oil 
% extraction, borrowing and default policies.
% Default decision as in Eaton and Gersovitz with exogenous reentry
%
% NOTE: you have to execute maduro_oil_run.m first
  
%% Simulation
   N  = 500;    % number of simulations
   T  = 500;    % Sample simulation Total number of periods (=T+J)
   J  = 1000;   % Drop first J simulations
   L  = 1;      % Dynamic lags for cross-correlations
   w  = 1600;   % H-P smoothing parameter
   nv = 7;      % Number of variables to compute sample moments
     
   M  = zeros(nv,1,N);      % Empty matrix to store the mean of the variables 
   SD = zeros(nv,1,N);      % Empty matrix to store the standard deviations of the variables
   XC = zeros(nv,1+2*L,N);  % Empty matrix to store the cross-corr of the variables
   AC = zeros(1+2*L,nv,N);  % Empty matrix to store the autocorrelation of the variables

   for i = 1:N

       st_init = simulmarkov(Pstar,T+J,2);
       
       st     = st_init(J+1:end); 
       xt     = x(st);             
       bp     = b(xt);            
       Yt     = Y(st);             
       yt     = log(Yt);
       Bt     = B(st);             
       qt     = min(max(0,funeval(cq,qspace,bp)),(1/(1+rate)));   
       Ct     = Yt+Bt-qt.*bp;               
       ct     = log(Ct);
       TBt    = Yt - Ct;                    
       TBt2Y  = TBt./Yt;
       
       bp_bar     = hpfilter(bp,w);
       Yt_bar     = hpfilter(Yt,w);
       yt_bar     = hpfilter(ct,w);
       Bt_bar     = hpfilter(Bt,w);
       Ct_bar     = hpfilter(Ct,w);
       ct_bar     = hpfilter(ct,w);
       TBt_bar    = hpfilter(TBt,w);
       TBt2Y_bar  = hpfilter(TBt2Y,w);
       qt_bar     = hpfilter(qt,1);
       
       corrcy         = corrcoef(ct_bar,yt_bar);
       xccy(:,:,i)    = corrcy; 
       corrTB2Yy      = corrcoef(TBt2Y_bar,yt_bar);
       xcTB2Yy(:,:,i) = corrTB2Yy; 
       [M(:,:,i),SD(:,:,i),XC(:,:,i),AC(:,:,i)]=samplemoms([Yt_bar yt_bar Ct_bar ct_bar TBt_bar TBt2Y_bar qt_bar],1,L);
      
   end

 % Take averages of the N-simulations over the third dimensions
   
   mn    = mean(M,3) ;
   sd    = mean(SD,3);
   xcorr = mean(XC,3);
   acorr = mean(AC,3);
   xcorrtb = mean(xcTB2Yy,3);
   xcorrcy = mean(xccy,3);
   % Plots the last simulation 
   
   figure;
   subplot(2,2,1)
   plot(Yt); title('Income','interpreter','latex');
   subplot(2,2,2)
   plot(Ct); title('Consumption','interpreter','latex');  
   subplot(2,2,3)
   plot(TBt2Y); title('Trade Balance to GDP','interpreter','latex');
   
   figure;
   subplot(1,3,1)
   spy(Pstar); title('Optimal Transition Matrix P(b,z)','interpreter','latex');
   subplot(1,3,2)
   plot(Eq); title('Equilibrium Bond Price','interpreter','latex');
   subplot(1,3,3)
   plot(pie); title('Ergodic distribution','interpreter','latex');

    figure; 
    subplot(2,2,1)
    plot(Yt_bar); title(' Detrended Income','interpreter','latex');
    subplot(2,2,2)
    plot(Ct_bar); title('Detrended Consumption','interpreter','latex');
    subplot(2,2,3)
    plot(TBt2Y_bar); title(' Detrended Trade Balance to GDP','interpreter','latex');
    
    
vp     = reshape(vpay,n2,n1);
vpd    = vp*mn_g;
vd     = reshape(vdef,n2,n1);
vdd    = vd*mn_g;
difv   = vp-vd;
diffvd = vpd-vdd;

if length(z)>1
    
ds     = reshape(D,n2,n1);
z2 =[min(z);z; max(z); max(z)]; % Useful when filling

 for i=1:n1
hgt    = (sum(ds(:,i)));
yll(i) = (((bmax-bmin)/n2)*hgt'+bmin); 
 end 
yll  = yll';
yll2 = [bmax; yll; yll(end); bmax];

figure('color', 'w');
subplot(1,2,1);
fill(yll2, z2, [0.9,0.9,0.9]);
xlim([bmin bmax]);
ylim([min(z) max(z)]);
set(gca,'Ytick',linspace(min(z),max(z),4),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Z','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

subplot(1,2,2);
plot(diffvd(52,:));
set(gca,'fontsize', 11);
legend('$V^{pay}$-$V^{def}$','Location','southoutside', 'Orientation', 'horizontal')
title(['Difference between Value Function under repayment and default'],'interpreter','latex');
xlabel('Z state','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

else 
  
    ds     = reshape(D,n2,n1);
g2 =[min(g);g; max(g); max(g)]; % Useful when filling

 for i=1:n1
hgt    = (sum(ds(:,i)));
yll(i) = (((bmax-bmin)/n2)*hgt'+bmin); 
 end 
yll  = yll';
yll2 = [bmax; yll; yll(end); bmax];

figure('color', 'w');
subplot(1,2,1);
fill(yll2, g2, [0.9,0.9,0.9]);
xlim([bmin bmax]);
ylim([min(g) max(g)]);
set(gca,'Ytick',linspace(min(g),max(g),4),'fontsize', 11);
hold on
fill(0,0,'w');
hold off
legend('Repayment set','Default set','Location','southoutside', 'Orientation', 'horizontal')
title(['Decision set under default and repayment'],'interpreter','latex');
xlabel('Assets $b_{t}$','interpreter','latex','FontSize',12)
ylabel('Z','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')

subplot(1,2,2);
plot(diffvd(52,:));
set(gca,'fontsize', 11);
legend('$V^{pay}$-$V^{def}$','Location','southoutside', 'Orientation', 'horizontal')
title(['Difference between Value Function under repayment and default'],'interpreter','latex');
xlabel('Z state','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
legend('boxoff')
end

   drs = ((1./qt-1)-(rate))*100;

   time = toc;   
%% Print summary
   
fprintf('\n')
fprintf('\n                    Moments             \n')
fprintf('                                  Model  ')
fprintf('\nVariable               Std.Dev Corr(i,p) Autocorr')
fprintf('\nIncome               %8.2f %8.2f %8.2f',sd(1),xcorr(1,L+1),acorr(L,1))
fprintf('\nLog. Income          %8.2f %8.2f %8.2f',sd(2),xcorr(2,L+1),acorr(L,2))
fprintf('\nConsumption          %8.2f %8.2f %8.2f',sd(3),xcorr(3,L+1),acorr(L,3))
fprintf('\nLog. Consumption     %8.2f %8.2f %8.2f',sd(4),xcorr(4,L+1),acorr(L,4))
fprintf('\nTrade Balance        %8.2f %8.2f %8.2f',sd(5),xcorr(5,L+1),acorr(L,5))
fprintf('\nTrade Balance to GDP %8.2f %8.2f %8.2f',sd(6),xcorr(6,L+1),acorr(L,6))
fprintf('\nBond Price           %8.2f %8.2f %8.2f',sd(7),xcorr(7,L+1),acorr(L,7))
fprintf('\nElapsed time (sec)   %8.2f\n',time)

