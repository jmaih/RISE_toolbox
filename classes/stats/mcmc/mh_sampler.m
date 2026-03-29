%  MH_SAMPLER -- Metropolis Hastings sampler
% 
%  ::
% 
% 
%    [Results]=MH_SAMPLER(logf,lb,ub)
% 
%    [Results]=MH_SAMPLER(logf,lb,ub,options)
% 
%    [Results]=MH_SAMPLER(logf,lb,ub,options,x0)
% 
%    [Results]=MH_SAMPLER(logf,lb,ub,options,x0,SIG)
% 
%  Args:
% 
%     - **logf** [char|function_handle]: Objective function to MINIMIZE!!!
% 
%     - **lb** [d x 1 vector]: lower bound of the paramters
% 
%     - **ub** [d x 1 vector]: upper bound of the paramters
% 
%     - **options** [struct]: see help for rsamplers.rwmh.rwmh
% 
%     - **x0** [d x 1 vector|empty|{.5*(lb+ub)}]: initial condition for the
%     sampler 
% 
%     - **SIG** [d x d|empty|{eye(d)}]: initial covariance matrix
% 
%  Returns:
% 
%     - **Results** [struct]:
%       - **pop** [nchain x N struct]: fields are "x" for the parameter vector
%       and "f" for the value of the parameter vector
%       - **best** [nchain x 1]: vector of best individual in each chain
%       - **m** [vector]: mean of the parameter draws
%       - **SIG** [matrix]: covariance of the parameter draws
%       - **c** [scalar]: final value of the scale on the covariance matrix
%       - **funevals** [1 x nchain vector]: function evaluations
%       - **stats** [struct]: stats on the optimization
% 
%  Note:
% 
%     - It is assumed that logf is a function to minimize
% 
%  Example:
% 
%     See also:	CONSTANT_BVAR_SAMPLER, RRF_SAMPLER
%