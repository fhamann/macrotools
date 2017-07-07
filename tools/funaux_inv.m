function EE = funaux_inv(func,ss,e,sp,z,w,beta,params)

% FUNAUX_INV  Euler Equation Function to be solved in a sparse
% representation. Uses the model function in which the inverse of fx and 
% g are specified. See: SOLVETI
% 
% Usage:
%                EE = funaux_inv(func,ss,e,sp,z,w,beta,params)
% INPUTS
%  func    : model function (see solveti)
%  ss      : endogenous state variable
%  e       : exogenous state variable
%  sp      : next period endogenous state variable
%  z       : initial policy function guess in a sparse representation
%  w       : transition probability matrix
%  beta    : model discount factor
%  params  : model parameters
%
% OUTPUTS
% EE       : Euler Equation 

  ss    = max(ss,min(z.range)+eps);
  sp    = max(sp,min(z.range)+eps);
  spvec = sp*ones(length(w),1);
   
  [~,xx]       = feval(func,'i',ss,[],e,sp,sp,params{:}); xx = max(xx,eps);
  [~,xxp]      = feval(func,'i',sp,[],e,spinterp(z,spvec),spinterp(z,spvec),params{:}); xxp = max(xxp,eps);
  [f,fx,fs]    = feval(func,'f',ss,xx,e,[],[],params{:}); 
  [fp,fxp,~]   = feval(func,'f',sp,xxp,e,[],[],params{:});
  [~,~,fsp]    = feval(func,'f',spinterp(z,sp),xxp,e,[],[],params{:});
  [g,gx,gs]    = feval(func,'g',sp,xx,e,[],[],params{:});
  [gp,gxp,gsp] = feval(func,'g',sp,xxp,e,[],[],params{:});
  
  EE = -fx./gx - (beta*w*(fsp-gs.*(fxp./gxp)));

end