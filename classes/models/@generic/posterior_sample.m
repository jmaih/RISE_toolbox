%--- help for generic/posterior_sample ---
%
%  Computes a sample of any quantity of interest using
%  parameter draws from a population e.g. a posterior simulation
% 
%  ::
% 
%    [result]=posterior_sample(m,pop,dowhat)
% 
%    [result]=posterior_sample(m,pop,dowhat,howmany)
% 
%    [result]=posterior_sample(m,pop,dowhat,howmany,ouf)
% 
%    [result]=posterior_sample(m,pop,dowhat,howmany,ouf,varargin)
% 
%    [result,is_failed,time_it_took]=posterior_sample(...)
% 
%  Args:
% 
%     - **m** [rise|dsge|svar|rfvar|valid rise object]: model object
% 
%     - **pop** [m x n struct]: parameter draws, with "m" the number of chains
%       and "n" the number of draws in each chain. Each element of "pop" is a
%       structure with fields "f" (not used), the value of the posterior and "x"
%       the parameter vector
% 
%     - **dowhat** [fhandle]: function (handle) to apply to each parameterized
%       model object. e.g. dowhat=@irf, dowhat=@simulate, dowhat=@forecast, etc.
%       "dowhat" need not be a method of "m": it represents the quantity of
%       interest.
% 
%     - **howmany** [integer|{m x n}]: number of draws to use in the
%       calculation
% 
%     - **ouf** [fhandle|{[]}]: output update function. Function that updates
%       the output before storing it. e.g. if dowhat=@filter, one may be
%       interested in the filters only and in that case ouf=@(x)x.filtering.
% 
%     - **varargin** [pairwise args]: valid pairwise arguments for the model
%       object
% 
%  Returns:
%     :
% 
%     - **result** [1 x howmany cell]: container of the various applications of
%       the "dowhat" handle
% 
%     - **time_it_took** [numeric]: number of seconds needed to run all the
%       simulations.
% 
%  Note:
% 
%     - the function will exploit parallel computation if there are workers
%       idle.
% 
%     - Because the solving of the model is sometimes iterative, a change of
%       solver or of the settings of the solver can result in the model not
%       being solved or more generally simulations failures. The algorithm will
%       loop until the requested number of simulations is obtained. But it will
%       not point to the parameter vectors that fail.
% 
%