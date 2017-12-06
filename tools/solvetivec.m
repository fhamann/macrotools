% SOLVETI   Solves stochastic dynamic program by time-iteration
% USAGE:
%                [x,sp,v] = solveti(model,fspace,s,xex)
%
% INPUTS
%  model   : structured array with model specifications (see CompEcon)
%  fspace  : function space (example: fspace = fundefn('lin',n,smin,smax)
%  s       : ns by 1 vector of state variable nodes
%  xex     : exogenous variable
%
% OUTPUTS
%  x       : ns by ne optimal controls
%  sp      : ns by ne optimal next period state
%  v       : ns by ne value function
%  resid   : nplot by ne residuals
%  Fopt    : ns by ne optimal reward function
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

% function [x,sp,v,resid,Fopt] = solveti(model,fspace,s,xex,b,y)

function [x,sp,Fopt,resid] = solvetivec(model,fspace,s,xex,b,y)

% SET THE DEFAULT OPTIONS
  tol_it  = tiget('solveti','tol_it', 10^-8);    % error tolerance for time iteration
  lambda  = tiget('solveti','lambda',0.8);       % updating weight for time iteration
  maxit   = tiget('solveti','maxit',10000);      % maximum iters.  for time iteration
  nplot   = tiget('solveti','nplot',000);       % # of nodes to evaluate errors (0 = do not do errors)
  prntit  = tiget('solveti','prntit',0);        % prints every prntit iters. (0 = do not print)

  if ~isfield(model,'e'); model.e = 0; e = 0;       else e = model.e; end;
  if ~isfield(model,'w'); model.w = eye(length(e)); else w = model.w;  end
  if ~isfield(model,'discount'); warning('model must have a discount factor');
  else beta = model.discount; end

func   = model.func;
params = model.params;

if nargin < 4
 % DETERMINE NUMBER OF DIMENSIONS & COORDINATES
  ns   = fspace.n;                          % Number of endogenous state variable nodes 
  ne   = length(e);                         % Number of exogenous state variable nodes 
  n    = ns*ne;                             % Number of state nodes
  npar = length(params);                    % Number of parameters
  smin = fspace.a;                          % Lower endogenous state variable limit
  smax = fspace.b;                          % Higher endogenous state variable limit
    
  if ~isfield(model,'actions'); sp = s; else sp = model.actions; end
  if nplot == 0 && nargout>3 ; display(['Error analysis can not be computed, set the number of error-evaluation nodes higher than zero']);
  return; end; 

% MATRIX REPLICATION
 ss = repmat(s,1,ne);                       % replicate the endogenous state variable 
 ee = repmat(e,ns,1);                       % replicate the exogenous state variable
 ww = repmat(w',length(s)/length(w),1);                      % replicate the transition probability matrix 
     
 snodes = funnode(fspace);               % state collocaton nodes

% TIME ITERATION ALGORITHM  
   if prntit>0;    disp('Iteration   Norm'); end
     
% COMPUTE BOUNDS
   [xl,xu] = feval(func,'b',s,[],e,sp,sp,params{:});

   for it = 1:maxit                                                          % perform iterations
    spold = sp;                                                              % store old policy function values
    [~,x]    = feval(func,'i',s,[],e,sp,sp,params{:});                    % find the control variable values given the stored policy function
    x = min((max(x,xl)),xu); 
    
    [f,fx,fs] = feval(func,'f',s,x,e,[],[],params{:});                    % evaluate the reward function 
    [g,gx,gs] = feval(func,'g',s,x,e,[],[],params{:});                    % evaluate the transition function    
    mu  = -fx./gx;                                                           % calculate the Lagrange multiplier using the reward and the transition function
    mup = minterp(snodes,mu,sp);                                % interpolate the Lagrange multiplier
    GSP = minterp(snodes,gs,sp);                % interpolate the transition function derivative according to the state given the policy function
    FSP = minterp(snodes,fs,sp);                % interpolate the reward function derivative according to the state given the policy function
    aux = (FSP + (GSP.*mup));                                               % create an auxiliar variable that incorporates the previous interpolation results
    aux = reshape(aux,length(s)/length(w),length(w));
    emu = reshape((w*aux')',length(s),1);                  % calculate the expected value of the Lagrange multiplier
    [xopt,~]   = feval(func,'i',s,x,e,-beta.*gx.*emu,sp,params{:});        % find the optimal control variable values
    xopt = min((max(xopt,xl)),xu);           
    [sp,gx,gs] = feval(func,'g',s,xopt,e,[],[],params{:});                 % find the policy function values given the optimal control variable values
     sp = min(max(sp,smin),smax);
    error_it   = max(max(abs(sp-spold)))/max(max(abs(spold)));               % compute the iteration errors

    if mod(it, prntit) == 0; fprintf('%5i   %10.1e\n',it,error_it); end % print the iteration error each prntit iterations (see tiset)
    if abs(error_it)<tol_it; break; end                                      % convergence check
 
    sp = lambda*sp + (1-lambda)*spold;                                       % update the policy function
   end
  if prntit>0; fprintf('%5i   %10.1e\n',it,error_it); end                   % Print last iteration result
  
  EE    = mu - beta.*emu;                                                     % compute the Euler Equation values
  [~,x] = feval(func,'i',s,[],e,sp,sp,params{:}); x = min((max(x,xl)),xu); % find the optimal control variable values 
  Fopt     = feval(func,'f',s,x,e,[],[],params{:});                           % evaluate the reward function
%   v     = Fopt/(eye(length(w))-beta*w);                                                % calculate the optimal value function
  
  % ERROR ANALYSIS
  if nplot > 0    
      
   splot  = linspace(smin,smax,nplot)';                                      % create an equidistant grid to evaluate the errors
   eplot  = linspace(min(e),max(e),nplot)';
   smat   = reshape(s,ns,length(w));
   emat   = reshape(s,ns,length(w));
   SS     = minterp(snodes,smat(:,1),splot);                              % endogenous state variable matrix replication
   E      = minterp(snodes,emat(:,1),eplot);                              % exogenous state variable matrix replication
   SP     = minterp(snodes,sp,splot);                                      % interpolate the error evaluation grid in the optimal policy function values
   if length(beta)>1
   BETA   = minterp(snodes,beta,splot);
   else 
   BETA   = beta;
   end
   [XL,XU] = feval(func,'b',SS,[],E,SP,SP,params{:});
   
   [~,XX] = feval(func,'i',SS,[],E,SP,SP,params{:});
   XX     = min((max(XX,XL)),XU);
   
   [F,Fx,Fs]  = feval(func,'f',SS,XX,E,[],[],params{:});
   [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,[],[],params{:});
  
   MU    = -Fx./Gx;
   MUP   = minterp(splot,MU,Sp);
   GSP   = minterp(splot,Gs,Sp);
   FSP   = minterp(splot,Fs,Sp);
   AUX   = repmat((FSP + (MUP.*GSP)),1,length(w));
   EMU   = (sum(w*AUX'))';
   resid = MU - BETA.*EMU;                                                    % compute the residual
  end
%%     
elseif nargin < 5
% DETERMINE NUMBER OF DIMENSIONS & COORDINATES
  ns   = fspace.n;                          % Number of endogenous state variable nodes 
  ne   = length(e);                         % Number of exogenous state variable nodes 
  n    = ns*ne;                             % Number of state nodes
  npar = length(params);                    % Number of parameters
  smin = fspace.a;                          % Lower endogenous state variable limit
  smax = fspace.b;                          % Higher endogenous state variable limit
    
  if ~isfield(model,'actions'); sp = repmat(s,1,ne); else sp = model.actions; end
  if nplot == 0 && nargout>3 ; display(['Error analysis can not be computed, set the number of error-evaluation nodes higher than zero']);
  return; end; 

% MATRIX REPLICATION
%  ss = repmat(s,1,ne);                       % replicate the endogenous state variable 
%  ee = repmat(e,ns,1);                       % replicate the exogenous state variable
%  ww = repmat(w',ns,1);                      % replicate the transition probability matrix 
    
% TIME ITERATION ALGORITHM  
    if prntit>0;    disp('Iteration   Norm'); end
     
% COMPUTE BOUNDS
   [xl,xu] = feval(func,'b',s,[],e,sp,sp,xex,params{:});

   for it = 1:maxit                                                          % perform iterations
    spold = sp;                                                              % store old policy function values
    [~,x]    = feval(func,'i',s,[],e,sp,sp,xex,params{:});                    % find the control variable values given the stored policy function
    x = min((max(x,xl)),xu); 
    
    [f,fx,fs] = feval(func,'f',s,x,e,[],[],xex,params{:});                    % evaluate the reward function 
    [g,gx,gs] = feval(func,'g',s,x,e,[],[],xex,params{:});                    % evaluate the transition function    
   
     mu  = -fx./gx;                                                           % calculate the Lagrange multiplier using the reward and the transition function
    mup = minterp(snodes,mu,sp);                                % interpolate the Lagrange multiplier
    GSP = minterp(snodes,gs,sp);                % interpolate the transition function derivative according to the state given the policy function
    FSP = minterp(snodes,fs,sp);                % interpolate the reward function derivative according to the state given the policy function
    aux = (FSP + (GSP.*mup));                                               % create an auxiliar variable that incorporates the previous interpolation results
    aux = reshape(aux,length(s)/length(w),length(w));
    emu = reshape((w*aux')',length(s),1);                  % calculate the expected value of the Lagrange multiplier
    [xopt,~]   = feval(func,'i',s,x,e,-beta*gx.*emu,sp,xex,params{:});        % find the optimal control variable values
    xopt = min((max(xopt,xl)),xu);           
    [sp,gx,gs] = feval(func,'g',s,xopt,e,[],[],xex,params{:});                 % find the policy function values given the optimal control variable values
     sp = min(max(sp,smin),smax);
    error_it   = max(max(abs(sp-spold)))/max(max(abs(spold)));                   % compute the iteration errors

    if mod(it, prntit) == 0; fprintf('%5i   %10.1e\n',it,error_it); end % print the iteration error each prntit iterations (see tiset)
    if abs(error_it)<tol_it; break; end                                      % convergence check
 
    sp = lambda*sp + (1-lambda)*spold;                                       % update the policy function
   end
  if prntit>0; fprintf('%5i   %10.1e\n',it,error_it); end                                % Print last iteration result
  EE    = mu - beta*emu;                                                     % compute the Euler Equation values
  [~,x] = feval(func,'i',s,[],e,sp,sp,xex,params{:}); x = min((max(x,xl)),xu); % find the optimal control variable values 
  Fopt  = feval(func,'f',s,x,e,[],[],xex,params{:});                           % evaluate the reward function
%   v     = Fopt/(eye(ne)-beta*w);                                                % calculate the optimal value function
  
  % ERROR ANALYSIS
  if nplot > 0    
      
   splot  = linspace(smin,smax,nplot)';                                      % create an equidistant grid to evaluate the errors
   SS     = minterp(snodes,s,splot);                                            % endogenous state variable matrix replication
   E      = minterp(snodes,e,splot);                                                % exogenous state variable matrix replication
   SP     = minterp(snodes,sp,splot);                           % interpolate the error evaluation grid in the optimal policy function values
   XEX    = minterp(snodes,xex,splot);
   if length(beta)>1
   BETA   = minterp(snodes,beta,splot);
   else 
   BETA   = beta;
   end
   [XL,XU] = feval(func,'b',SS,[],E,SP,SP,XEX,params{:});
   
   [~,XX] = feval(func,'i',SS,[],E,SP,SP,XEX,params{:});
   XX     = min((max(XX,XL)),XU);
   
   [F,Fx,Fs]  = feval(func,'f',SS,XX,E,[],[],XEX,params{:});
   [Sp,Gx,Gs] = feval(func,'g',SS,XX,E,[],[],XEX,params{:});
  
   MU    = -Fx./Gx;
   MUP   = minterp(splot,MU,Sp);
   GSP   = minterp(splot,Gs,Sp);
   FSP   = minterp(splot,Fs,Sp);
   AUX   = repmat((FSP + (MUP.*GSP)),1,length(w));
   EMU   = (sum(w*AUX'))';
   resid = MU - BETA.*EMU;                                                    % compute the residual   
   end
end
    
