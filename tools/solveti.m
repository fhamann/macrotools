% SOLVETI   Solves stochastic dynamic program by time-iteration
% USAGE:
%                [x,sp,v] = solveti(model,fspace,s,xex)
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
% MODEL STRUCTURE FIELDS
%   func           : function file name
%   discount       : discount factor
%   e              : exogenous state variable
%   w              : transition probability matrix    
%   actions        : vector of discrete actions
%   params         : additional parameters to function file
% FUNCTION FILE FORMAT
%   [out1,out2,out3] = func(flag,s,x,e,additional parameters)
%   if flag = 'f' returns reward function and derivatives
%      f, fx, fs
%   if flag = 'g' returns transition function and derivatives
%      g, gx, gs
%   if flag = 'i' returns inverse functions of fx and g
%      fxi, gi
%
% by Franz Hamann and Juan Camilo Mendez based on Fackler and Miranda CompEcon Toolbox.

function [x,sp,v,resid,s] = solveti(model,fspace,s,xex)

% SET THE DEFAULT OPTIONS
  tol_it  = tiget('solveti','tol_it', 10^-6);    % error tolerance for time iteration
  lambda  = tiget('solveti','lambda',0.5);       % updating weight for time iteration
  maxit   = tiget('solveti','maxit',10000);      % maximum iters.  for time iteration
  nplot   = tiget('solveti','nplot',2000);       % # of nodes to evaluate errors (0 = do not do errors)
  prntit  = tiget('solveti','prntit',25);        % prints every prntit iters. (0 = do not print)

  if ~isfield(model,'e'); model.e = 0; e = 0;       else e = model.e'; end;
  if ~isfield(model,'w'); model.w = eye(length(e)); else w = model.w;  end
  if ~isfield(model,'discount'); warning('model must have a discount factor');
  else beta = model.discount; end

func   = model.func;
params = model.params;

% DETERMINE NUMBER OF DIMENSIONS & COORDINATES
  ns   = fspace.n;                          % Number of endogenous state variable nodes 
  ne   = length(e);                         % Number of exogenous state variable nodes 
  n    = ns*ne;                             % Number of state nodes
  npar = length(params);                    % Number of parameters
  smin = fspace.a;                          % Lower endogenous state variable limit
  smax = fspace.b;                          % Higher endogenous state variable limit
  
  if ~isfield(model,'actions'); sp = repmat(s,1,ne); else sp = model.actions; end
  if size(sp,1)>ns; warning('actions must have the same length of states');   end
  if nplot == 0 && nargout>3 ; display(['Error analysis can not be computed, set the number of error-evaluation nodes higher than zero']);
  return; end; 

% MATRIX REPLICATION
 ss = repmat(s,1,ne);                       % replicate the endogenous state variable 
 ee = repmat(e,ns,1);                       % replicate the exogenous state variable
 ww = repmat(w',ns,1);                      % replicate the transition probability matrix 

if nargin < 4
% TIME ITERATION ALGORITHM  
<<<<<<< Updated upstream
   if (nargin(func)-npar) > 4                                                % MODEL WITH THE SPECIFICATION OF THE INVERSE OF FX AND G 
=======
   if (nargin(func)-npar) > 5                                                  % MODEL WITH THE SPECIFICATION OF THE INVERSE OF FX AND G 
>>>>>>> Stashed changes
   disp('Iteration   Norm');
     
% COMPUTE BOUNDS
   [xl,xu] = feval(func,'b',ss,[],ee,sp,sp,params{:});

   for it = 1:maxit                                                          % perform iterations
    spold = sp;                                                              % store old policy function values
    [~,xx]    = feval(func,'i',ss,[],ee,sp,sp,params{:});                    % find the control variable values given the stored policy function
    xx = min((max(xx,xl)),xu); 
    
    [f,fx,fs] = feval(func,'f',ss,xx,ee,[],[],params{:});                    % evaluate the reward function 
    [g,gx,gs] = feval(func,'g',ss,xx,ee,[],[],params{:});                    % evaluate the transition function    
   
    mu  = -fx./gx;                                                           % calculate the Lagrange multiplier using the reward and the transition function
    mup = interp1(s,mu,sp,'linear','extrap');                                % interpolate the Lagrange multiplier

    MUP = reshape(mup,n,ne);                                                 % make the Lagrange multiplier interpolant conformable
    GSP = reshape((interp1(s,gs,sp,'linear','extrap')),n,ne);                % interpolate the transition function derivative according to the state given the policy function
    FSP = reshape((interp1(s,fs,sp,'linear','extrap')),n,ne);                % interpolate the reward function derivative according to the state given the policy function
    aux = (FSP + (GSP.*MUP))';                                               % create an auxiliar variable that incorporates the previous interpolation results
    Emu = reshape(sum(reshape(ww(:).*aux(:),ne,n)),ns,ne);                   % calculate the expected value of the Lagrange multiplier
    [xopt,~]   = feval(func,'i',ss,xx,ee,-beta*gx.*Emu,sp,params{:});        % find the optimal control variable values
    xopt = min((max(xopt,xl)),xu);           
    [sp,gx,gs] = feval(func,'g',ss,xopt,ee,[],[],params{:});                 % find the policy function values given the optimal control variable values
     sp = min(max(sp,smin),smax);
    error_it   = max(max(abs(sp-spold)))/max(max(abs(spold)));               % compute the iteration errors

    if mod(it, prntit) == 0; fprintf('%5i   %10.1e\n',it,error_it); end % print the iteration error each prntit iterations (see tiset)
    if abs(error_it)<tol_it; break; end                                      % convergence check
 
    sp = lambda*sp + (1-lambda)*spold;                                       % update the policy function
   end
  fprintf('%5i   %10.1e\n',it,error_it);                                % Print last iteration result
  EE    = mu - beta*Emu;                                                     % compute the Euler Equation values
  [~,x] = feval(func,'i',ss,[],ee,sp,sp,params{:}); x = min((max(x,xl)),xu); % find the optimal control variable values 
  F     = feval(func,'f',ss,x,ee,[],[],params{:});                           % evaluate the reward function
  v     = F/(eye(ne)-beta*w);                                                % calculate the optimal value function
  
  % ERROR ANALYSIS
  if nplot > 0    
      
   splot  = linspace(smin,smax,nplot)';                                      % create an equidistant grid to evaluate the errors
   SS     = repmat(splot,1,ne);                                              % endogenous state variable matrix replication
   E      = repmat(e,nplot,1);                                               % exogenous state variable matrix replication
   W      = repmat(w',nplot,1);                                              % transition probability matrix replication
   SP     = interp1(s,sp,splot,'linear','extrap');                           % interpolate the error evaluation grid in the optimal policy function values
      
   [XL,XU] = feval(func,'b',SS,[],E,SP,SP,params{:});
   
   [~,XX] = feval(func,'i',SS,[],E,SP,SP,params{:});
   XX     = min((max(XX,XL)),XU);
   
   [F,Fx,Fs]  = feval(func,'f',SS,XX,E,[],[],params{:});
   [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,[],[],params{:});
  
   MU    = -Fx./Gx;
   MUP   = reshape((interp1(splot,MU,Sp,'linear','extrap')),ne*nplot,ne);
   GSP   = reshape((interp1(splot,Gs,Sp,'linear','extrap')),ne*nplot,ne);
   FSP   = reshape((interp1(splot,Fs,Sp,'linear','extrap')),ne*nplot,ne);
   AUX   = (FSP + (MUP.*GSP))';
   EMU   = reshape(sum(reshape(W(:).*AUX(:),ne,nplot*ne)),nplot,ne);    
   resid = MU - beta*EMU;                                                    % compute the residual
  end

 elseif (nargin(func)-npar) <= 5                                             % MODEL WITHOUT THE SPECIFICATION OF THE INVERSE OF FX AND G (uses Matlab symbolic toolbox)
%  COMPUTE BOUNDS
  [xl,xu] = feval(func,'b',ss,[],ee,params{:});   
     
  syms ssym xsym esym                                                        % convert into symbols the states variables and the control
  gsym = symfun(feval(func,'g',ssym,xsym,esym,params{:}),[ssym xsym esym]);  % express the transition function as a symbolic function 
  gi   = matlabFunction(finverse(gsym,xsym));                                % rewrite the transition function in terms of the policy function.

  [~,fxsym,~] = feval(func,'f',ssym,xsym,esym,params{:});                    % express the reward function and its derivatives as a symbolic function 
  fxi = matlabFunction(finverse(fxsym,xsym));                                % rewrite the reward function derivative in terms of the policy function.
  disp('Iter         Norm');
  
  for it = 1:maxit
   spold = sp;
   xx = feval(gi,ss,sp,ee);       
   xx = min((max(xx,xl)),xu);

   [f,fx,fs]  = feval(func,'f',ss,xx,ee,params{:});
   [sp,gx,gs] = feval(func,'g',ss,xx,ee,params{:});
   
   mu   = -fx./gx;
   mup  = interp1(s,mu,sp,'linear','extrap');
   MUP  = reshape(mup,n,ne);   
   GSP  = reshape((interp1(s,gs,sp,'linear','extrap')),n,ne);
   FSP  = reshape((interp1(s,fs,sp,'linear','extrap')),n,ne);
   aux  = (FSP + (GSP.*MUP))';
   Emu  = reshape(sum(reshape(ww(:).*aux(:),ne,n)),ns,ne);
   if nargin(fxi) >1;  
     xopt = feval(fxi,ee,ss,-beta*gx.*Emu);
   else 
     xopt = feval(fxi,-beta*gx.*Emu);
   end;
   xopt =  min((max(xopt,xl)),xu);
   [sp,gx,~] = feval(func,'g',ss,xopt,ee,params{:}); 
  
   error_it  = max(max(abs(sp-spold)))/max(max(abs(spold)));

   if mod(it, prntit) == 0; fprintf('%5i   %10.1e\n',it,error_it); end
   if abs(error_it)<tol_it; break; end
   sp = lambda*sp + (1-lambda)*spold;
  
  end
  fprintf('%5i     %10.1e\n',it,error_it);    % Print last iteration result
  EE  = mu - beta*Emu;
  x   = feval(gi,ss,sp,ee); x = min((max(x,xl)),xu);
  F   = feval(func,'f',ss,x,ee,params{:});
  v   = F/(eye(ne)-beta*w);

% ERROR ANALYSIS  
  
  if nplot > 0   
   splot = linspace(smin*0.2,smax,nplot)'; 
   SS    = repmat(splot,1,ne);
   E     = repmat(e,nplot,1);
   SP    = interp1(s,sp,splot,'linear','extrap');
   SP    = min(max(SP,smin),smax);
   W     = repmat(w',nplot,1);
   XX    = feval(gi,SS,SP,E);
   XX     = min((max(XX,xl)),xu);
   [XL,XU] = feval(func,'b',SS,[],E,params{:}); 

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
 
else 
% TIME ITERATION ALGORITHM  
   if (nargin(func)-npar) > 5                                                  % MODEL WITH THE SPECIFICATION OF THE INVERSE OF FX AND G 
   disp('Iteration   Norm');
     
% COMPUTE BOUNDS
   [xl,xu] = feval(func,'b',ss,[],ee,sp,sp,xex,params{:});

   for it = 1:maxit                                                          % perform iterations
    spold = sp;                                                              % store old policy function values
    [~,xx]    = feval(func,'i',ss,[],ee,sp,sp,xex,params{:});                    % find the control variable values given the stored policy function
    xx = min((max(xx,xl)),xu); 
    
    [f,fx,fs] = feval(func,'f',ss,xx,ee,[],[],xex,params{:});                    % evaluate the reward function 
    [g,gx,gs] = feval(func,'g',ss,xx,ee,[],[],xex,params{:});                    % evaluate the transition function    
   
    mu  = -fx./gx;                                                           % calculate the Lagrange multiplier using the reward and the transition function
    mup = interp1(s,mu,sp,'linear','extrap');                                % interpolate the Lagrange multiplier

    MUP = reshape(mup,n,ne);                                                 % make the Lagrange multiplier interpolant conformable
    GSP = reshape((interp1(s,gs,sp,'linear','extrap')),n,ne);                % interpolate the transition function derivative according to the state given the policy function
    FSP = reshape((interp1(s,fs,sp,'linear','extrap')),n,ne);                % interpolate the reward function derivative according to the state given the policy function
    aux = (FSP + (GSP.*MUP))';                                               % create an auxiliar variable that incorporates the previous interpolation results
    Emu = reshape(sum(reshape(ww(:).*aux(:),ne,n)),ns,ne);                   % calculate the expected value of the Lagrange multiplier
    [xopt,~]   = feval(func,'i',ss,xx,ee,-beta*gx.*Emu,sp,xex,params{:});        % find the optimal control variable values
    xopt = min((max(xopt,xl)),xu);           
    [sp,gx,gs] = feval(func,'g',ss,xopt,ee,[],[],xex,params{:});                 % find the policy function values given the optimal control variable values
     sp = min(max(sp,smin),smax);
    error_it   = max(max(abs(sp-spold)))/max(max(abs(spold)));               % compute the iteration errors

    if mod(it, prntit) == 0; fprintf('%5i   %10.1e\n',it,error_it); end % print the iteration error each prntit iterations (see tiset)
    if abs(error_it)<tol_it; break; end                                      % convergence check
 
    sp = lambda*sp + (1-lambda)*spold;                                       % update the policy function
   end
  fprintf('%5i   %10.1e\n',it,error_it);                                % Print last iteration result
  EE    = mu - beta*Emu;                                                     % compute the Euler Equation values
  [~,x] = feval(func,'i',ss,[],ee,sp,sp,xex,params{:}); x = min((max(x,xl)),xu); % find the optimal control variable values 
  F     = feval(func,'f',ss,x,ee,[],[],xex,params{:});                           % evaluate the reward function
  v     = F/(eye(ne)-beta*w);                                                % calculate the optimal value function
  
  % ERROR ANALYSIS
  if nplot > 0    
      
   splot  = linspace(smin,smax,nplot)';                                      % create an equidistant grid to evaluate the errors
   SS     = repmat(splot,1,ne);                                              % endogenous state variable matrix replication
   E      = repmat(e,nplot,1);                                               % exogenous state variable matrix replication
   W      = repmat(w',nplot,1);                                              % transition probability matrix replication
   SP     = interp1(s,sp,splot,'linear','extrap');                           % interpolate the error evaluation grid in the optimal policy function values
      
   [XL,XU] = feval(func,'b',SS,[],E,SP,SP,xex,params{:});
   
   [~,XX] = feval(func,'i',SS,[],E,SP,SP,xex,params{:});
   XX     = min((max(XX,XL)),XU);
   
   [F,Fx,Fs]  = feval(func,'f',SS,XX,E,[],[],xex,params{:});
   [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,[],[],xex,params{:});
  
   MU    = -Fx./Gx;
   MUP   = reshape((interp1(splot,MU,Sp,'linear','extrap')),ne*nplot,ne);
   GSP   = reshape((interp1(splot,Gs,Sp,'linear','extrap')),ne*nplot,ne);
   FSP   = reshape((interp1(splot,Fs,Sp,'linear','extrap')),ne*nplot,ne);
   AUX   = (FSP + (MUP.*GSP))';
   EMU   = reshape(sum(reshape(W(:).*AUX(:),ne,nplot*ne)),nplot,ne);    
   resid = MU - beta*EMU;                                                    % compute the residual
  end

 elseif (nargin(func)-npar) <= 5                                             % MODEL WITHOUT THE SPECIFICATION OF THE INVERSE OF FX AND G (uses Matlab symbolic toolbox)
%  COMPUTE BOUNDS
  [xl,xu] = feval(func,'b',ss,[],ee,params{:});   
     
  syms ssym xsym esym xexsym                                                      % convert into symbols the states variables and the control
  gsym = symfun(feval(func,'g',ssym,xsym,esym,xexsym,params{:}),[ssym xsym esym xexsym]);  % express the transition function as a symbolic function 
  gi   = matlabFunction(finverse(gsym,xsym));                                % rewrite the transition function in terms of the policy function.

  [~,fxsym,~] = feval(func,'f',ssym,xsym,esym,xexsym,params{:});                    % express the reward function and its derivatives as a symbolic function 
  fxi = matlabFunction(finverse(fxsym,xsym));                                % rewrite the reward function derivative in terms of the policy function.
  disp('Iter         Norm');
  
  for it = 1:maxit
   spold = sp;
   xx = feval(gi,ss,sp,ee,xex);       
   xx = min((max(xx,xl)),xu);

   [f,fx,fs]  = feval(func,'f',ss,xx,ee,xex,params{:});
   [sp,gx,gs] = feval(func,'g',ss,xx,ee,xex,params{:});
   
   mu   = -fx./gx;
   mup  = interp1(s,mu,sp,'linear','extrap');
   MUP  = reshape(mup,n,ne);   
   GSP  = reshape((interp1(s,gs,sp,'linear','extrap')),n,ne);
   FSP  = reshape((interp1(s,fs,sp,'linear','extrap')),n,ne);
   aux  = (FSP + (GSP.*MUP))';
   Emu  = reshape(sum(reshape(ww(:).*aux(:),ne,n)),ns,ne);
   if nargin(fxi) >1;  
     xopt = feval(fxi,ee,ss,-beta*gx.*Emu,xex);
   else 
     xopt = feval(fxi,-beta*gx.*Emu);
   end;
   xopt =  min((max(xopt,xl)),xu);
   [sp,gx,~] = feval(func,'g',ss,xopt,ee,xex,params{:}); 
  
   error_it  = max(max(abs(sp-spold)))/max(max(abs(spold)));

   if mod(it, prntit) == 0; fprintf('%5i   %10.1e\n',it,error_it); end
   if abs(error_it)<tol_it; break; end
   sp = lambda*sp + (1-lambda)*spold;
  
  end
  fprintf('%5i     %10.1e\n',it,error_it);    % Print last iteration result
  EE  = mu - beta*Emu;
  x   = feval(gi,ss,sp,ee); x = min((max(x,xl)),xu);
  F   = feval(func,'f',ss,x,ee,xex,params{:});
  v   = F/(eye(ne)-beta*w);

% ERROR ANALYSIS  
  
  if nplot > 0   
   splot = linspace(smin*0.2,smax,nplot)'; 
   SS    = repmat(splot,1,ne);
   E     = repmat(e,nplot,1);
   SP    = interp1(s,sp,splot,'linear','extrap');
   SP    = min(max(SP,smin),smax);
   W     = repmat(w',nplot,1);
   XX    = feval(gi,SS,SP,E,xex);
   XX     = min((max(XX,xl)),xu);
   [XL,XU] = feval(func,'b',SS,[],E,xex,params{:}); 

   [F,Fx,Fs]  =  feval(func,'f',SS,XX,E,xex,params{:});
   [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,xex,params{:});

   MU    = -Fx./Gx;
   MUP   = reshape((interp1(splot,MU,Sp,'linear','extrap')),ne*nplot,ne);
   GSP   = reshape((interp1(splot,Gs,Sp,'linear','extrap')),ne*nplot,ne);
   FSP   = reshape((interp1(splot,Fs,Sp,'linear','extrap')),ne*nplot,ne);
   AUX   = (FSP+(GSP.*MUP))';
   EMU   = reshape(sum(reshape(W(:).*AUX(:),ne,nplot*ne)),nplot,ne);
   resid = MU - beta*EMU;
  end
   end 
end    
    
    
