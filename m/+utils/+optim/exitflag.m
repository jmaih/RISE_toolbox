%  INTERNAL FUNCTION: recast the exit flag from fsolve and lsqnonlin
% 
%  ::
% 
%    flag=exitflag(flag,x1,fx1,TolFun)
% 
%  Args:
% 
%     - **flag** [scalar]: exit flag from fsolve or lsqnonlin
%     - **x1** [scalar|vector|matrix]: optimum of the function to set to zero
%       (f(x)=0)
%     - **fx1** [scalar]: value or norm of f(x1)
%     - **TolFun** [{sqrt(eps)}|scalar]: scalar or vector of model objects
% 
%  Returns:
%     :
% 
%     - **flag** [scalar]: recomputed exitflag
% 
%