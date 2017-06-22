function [x,sp,v,resid] = solveti(model,fspace,s)
% SOLVETI   Solves stochastic dynamic program by time-iteration
%
% Usage:
%                [x,sp,v] = solveti(model,fspace,s)
%
% INPUTS
%  model   : structured array with model specifications (see CompEcon)
%  fspace  : function space (example: fspace = fundefn('lin',n,smin,smax)
%  s       : ns by 1 vector of state variable nodes
%
% OUTPUTS
%  x       : ns by ne optimal controls
%  sp      : ns by ne optimal next period state
%  v       : ns by ne value function
%  resid   : nplot by ne residuals
%
% where ne is the number of discretized exogenous states, e. 
%
% by Juan Camilo Méndez based on Fackler and Miranda CompEcon Toolbox.

%  Set the options
algor   = tiget('solveti','algor','fsolve');   % Euler Equation solution method
tol_it  = tiget('solveti','tol_it', 10^-12);   % error tolerance for time iteration
lambda  = tiget('solveti','lambda',0.5);       % updating weight for time iteration
maxit   = tiget('solveti','maxit',1000000);    % maximum iters.  for time iteration
nplot   = tiget('solveti','nplot',2000);       % # of nodes to evaluate errors (0 = do not do errors)
prntit  = tiget('solveti','prntit',20);        % prints every prntit iters. (0 = do not print)


if ~isfield(model,'e'); model.e = 0; e = 0;
    else e = model.e'; end;
if ~isfield(model,'w'); model.w = 1;
    else w = model.w; end
if ~isfield(model,'discount'); warning('model must have a discount factor');
    else beta = model.discount; end

func   = model.func;
params = model.params;

% Determine number of dimensions and coordinates
ns   = fspace.n;
ne   = length(e);
n    = ns*ne;
npar = length(params);
smin = fspace.a;
smax = fspace.b;

if ~isfield(model,'actions'); sp = repmat(0.9*smin+0.1*(smax+smin)/2,ns,ne);
    else sp = model.actions; end

if size(sp,1)>ns; warning('actions must have the same length of states');
end

if nplot == 0 && nargout>3 ; display(['Error analysis can not be computed, set the number of error-evaluation nodes higher than zero']);
return; end; 

ss = repmat(s,1,ne);
ee = repmat(e,ns,1);
ww = repmat(w',ns,1);
sp = ss;

if (nargin(func)-npar) > 4

for it = 1:maxit

  spold = sp;
  [~,xx]    = feval(func,'i',ss,[],ee,sp,sp,params{:}); xx = max(xx,0);
  

  [f,fx,fs] = feval(func,'f',ss,xx,ee,[],[],params{:});
  [g,gx,gs] = feval(func,'g',ss,xx,ee,[],[],params{:});

  
  mu  = -fx./gx;
  mup = interp1(s,mu,sp,'linear','extrap');

  MUP = reshape(mup,n,ne);   
  GSP = reshape((interp1(s,gs,sp,'linear','extrap')),n,ne);
  FSP = reshape((interp1(s,fs,sp,'linear','extrap')),n,ne);
  aux = (FSP + (GSP.*MUP))';


  Emu = reshape(sum(reshape(ww(:).*aux(:),ne,n)),ns,ne);

  [xopt,~] = feval(func,'i',ss,xx,ee,-beta*gx.*Emu,sp,params{:});
   xopt = max(0,xopt);

  [sp,gx,gs] = feval(func,'g',ss,xopt,ee,[],[],params{:}); 
 
  sp = max(sp,smin);
  
  error_it = max(max(abs(sp-spold)))/max(max(abs(spold)));

  if mod(it, prntit) == 0; fprintf('%d  %1.7f \n',it,error_it); end
  if abs(error_it)<tol_it; break; end
 
  sp = lambda*sp + (1-lambda)*spold;
end
   EE    = mu - beta*Emu;
   [~,x] = feval(func,'i',ss,[],ee,sp,sp,params{:}); x = max(x,0);
    u    = feval(func,'f',ss,x,ee,[],[],params{:});
    v    = u/(eye(ne)-beta*w);

if nplot > 0     % Error Analysis
    splot  = linspace(smin*0.25,smax,nplot)'; 
    SS     = repmat(splot,1,ne);
    E      = repmat(e,nplot,1);
    SP     = interp1(s,sp,splot);
    W      = repmat(w',nplot,1);
    [~,XX] = feval(func,'i',SS,[],E,SP,SP,params{:});
    XX     = max(XX,0);

    [F,Fx,Fs]  = feval(func,'f',SS,XX,E,[],[],params{:});
    [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,[],[],params{:});

    MU    = -Fx./Gx;
    MUP   = reshape((interp1(splot,MU,Sp,'linear','extrap')),ne*nplot,ne);
    GSP   = reshape((interp1(splot,Gs,Sp,'linear','extrap')),ne*nplot,ne);
    FSP   = reshape((interp1(splot,Fs,Sp,'linear','extrap')),ne*nplot,ne);
    AUX   = (FSP + (MUP.*GSP))';
    EMU   = reshape(sum(reshape(W(:).*AUX(:),ne,nplot*ne)),nplot,ne);
    resid = MU - beta*EMU;
 end

elseif (nargin(func)-npar) <= 4

syms ssym xsym esym
gsym = symfun(feval(func,'g',ssym,xsym,esym,params{:}),[ssym xsym esym]);
gi   = matlabFunction(finverse(gsym,xsym));

[fsym,fxsym,fssym] = feval(func,'f',ssym,xsym,esym,params{:});
fxi = matlabFunction(finverse(fxsym,xsym));

for it = 1:maxit

  spold = sp;

  xx = max(feval(gi,ss,sp,ee),0);

  [f,fx,fs]  = feval(func,'f',ss,xx,ee,params{:});
  [sp,gx,gs] = feval(func,'g',ss,xx,ee,params{:});
  mu  = -fx./gx;
  mup = interp1(s,mu,sp,'linear','extrap');
  MUP = reshape(mup,n,ne);   
  GSP = reshape((interp1(s,gs,sp,'linear','extrap')),n,ne);
  FSP = reshape((interp1(s,fs,sp,'linear','extrap')),n,ne);
  aux = (FSP + (GSP.*MUP))';
  Emu = reshape(sum(reshape(ww(:).*aux(:),ne,n)),ns,ne);

  xopt = feval(fxi,ee,ss,-beta*gx.*Emu);
  
  [sp,gx,gs] = feval(func,'g',ss,xopt,ee,params{:});
  sp = max(sp,smin);

  error_it = max(max(abs(sp-spold)))/max(max(abs(spold)));

  if mod(it, prntit) == 0; fprintf('%d  %1.7f \n',it,error_it); end
  if abs(error_it)<tol_it; break; end
 
  sp = lambda*sp + (1-lambda)*spold;
  
end
    EE  = mu - beta*Emu;
    x   = max(feval(gi,ss,sp,ee),0);
    u   = feval(func,'f',ss,x,ee,params{:});
    v   = u/(eye(ne)-beta*w);

if nplot > 0   % Error analysis
    splot = linspace(smin*0.2,smax,nplot)'; 
    SS    = repmat(splot,1,ne);
    E     = repmat(e,nplot,1);
    SP    = interp1(s,sp,splot,'linear','extrap');
    SP    = min(max(SP,smin),smax);
    W     = repmat(w',nplot,1);
    XX    = max(feval(gi,SS,SP,E),0);

    [F,Fx,Fs]  =  feval(func,'f',SS,XX,E,params{:});
    [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,params{:});

    MU    = -Fx./Gx;
    MUP   = reshape((interp1(splot,MU,Sp,'linear','extrap')),ne*nplot,ne);
    GSP   = reshape((interp1(splot,Gs,Sp,'linear','extrap')),ne*nplot,ne);
    FSP   = reshape((interp1(splot,Fs,Sp,'linear','extrap')),ne*nplot,ne);
    AUX   = (FSP+(GSP.*MUP))';
    EMU   = reshape(sum(reshape(W(:).*AUX(:),ne,nplot*ne)),nplot,ne);
    resid = MU - beta*EMU;
 end
end 



















