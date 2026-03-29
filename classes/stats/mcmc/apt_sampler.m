%  APT_SAMPLER -- Adaptive parallel tempering sampler
% 
%  ::
% 
% 
%    [Results]=APT_SAMPLER(logf,lb,ub)
% 
%    [Results]=APT_SAMPLER(logf,lb,ub,options)
% 
%    [Results]=APT_SAMPLER(logf,lb,ub,options,x0)
% 
%    [Results]=APT_SAMPLER(logf,lb,ub,options,x0,SIG)
% 
%  Args:
% 
%     - **logf** [char|function_handle]: Objective function to MINIMIZE!!!
% 
%     - **lb** [d x 1 vector]: lower bound of the paramters
% 
%     - **ub** [d x 1 vector]: upper bound of the paramters
% 
%     - **options** [struct]:
%       - **burnin** [integer|{0}]: number of burn-in initial simulations
%       - **N** [integer|{20000}]: number of simulations
%       - **c** [scalar|{0.25}]: initial scale for the covariance matrix
%       - **thinning** [integer|{1}]: number of thinning simulations. 1 means we
%       keep every draw, 2 means we keep every second, 3 every third, etc.
%       - **nchain** [integer|{1}]: number of parallel chains
%       - **MaxTime** [numeric|{inf}]: maximum simulation time.
%       - **MaxFunevals** [numeric|{inf}]: maximum number of function evaluations
%       - **proposal** [empty|'t-student'|{'normal'}] : proposal distribution.
%       - **alpha** [empty|{0.33}] : target acceptance rate
%       - **H** [empty|{12}] : number of striations or tempering stages
%       - **alpha_swap** [empty|{0.234}] : acceptance rate for jumps between
%       striations or target acceptance rate for swaps
%       - **rwm_fixed_p** [empty|{0}] : prob of drawing from a fixed initial
%       proposal 
%       - **rwm_exp** [empty|{0.6}] : Exponent of random-walk adaptation
%       step size
%       - **sw_exp_rm** [empty|{0.6}] : Exponent of temperature adaptation
%       step size 
%       - **ram_adapt** [true|{false}] : Use the robust AM adaptation 
%       - **separate_shape_adaptation** [false|{true}] : separate covariance
%       for each temperature 
%       - **fixed_temperatures** [true|{false}] : option for setting fixed
%       temperatures or not
% 
%     - **x0** [d x 1 vector|empty|{.5*(lb+ub)}]: initial condition for the
%     sampler 
% 
%     - **SIG** [d x d|empty|{eye(d)}]: initial covariance matrix
% 
%  Returns:
%     :
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