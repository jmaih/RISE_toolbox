function obj=estimate(obj,varargin)
% estimate - estimates the parameters of a RISE model
%
% Syntax
% -------
% ::
%
%   obj=estimate(obj)
%   obj=estimate(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **varargin** additional optional inputs among which the most relevant
%   for estimation are:
%
% - **estim_parallel** [integer|{1}]: Number of starting values
%
% - **estim_start_from_mode** [true|false|{[]}]: when empty, the user is
%   prompted to answer the question as to whether to start estimation from
%   a previously found mode or not. If true or false, no question is asked.
%
% - **estim_start_date** [numeric|char|serial date]: date of the first
%   observation to use in the dataset provided for estimation
%
% - **estim_end_date** [numeric|char|serial date]: date of the last
%   observation to use in the dataset provided for estimation
%
% - **estim_max_trials** [integer|{500}]: When the initial value of the
%   log-likelihood is too low, RISE uniformly draws from the prior support
%   in search for a better starting point. It will try this for a maximum
%   number of **estim_max_trials** times before squeaking with an error.
%
% - **estim_start_vals** [{[]}|struct]: when not empty, the parameters
%   whose names are fields of the structure will see their start values
%   updated or overriden by the information in **estim_start_vals**. There
%   is no need to provide values to update the start values for the
%   estimated parameters.
%
% - **estim_general_restrictions** [{[]}|function handle|cell array]: when
%   not empty, the argument can be a function handle or a cell array
%   containing the function handle and additional input arguments. The
%   general syntax for the calling the function handle is
%   viol=myfunc(obj,varargin), with **obj** the parameterized RISE object
%   which will be used in the computation of the restrictions violations.
%   Hence the restrictions are entered either as  @myfunc or as
%   {@myfunc,arg2,arg3,...}. Originally, RISE will call the function
%   without any inputs. In that case, RISE expects the output to be a
%   structure with fields :
%       - **number_of_restrictions** : number of restrictions
%       - **kf_filtering_level** [0|1|2|3]: if 0, no filters are required
%       for the computation of the restrictions. If 1, only the filtered
%       variables are required. If 2, the updated variables are required.
%       If 3, the smoothed variables are required.
%   When the function is called with inputs, RISE expects as output the
%   values of the restrictions. The sign of the violations does not matter.
%   All the user has to do is to put a zero where the restrictions are not
%   violated.
%
% - **estim_linear_restrictions** [{[]}|cell]: This is most often used in
%   the estimation of rfvar or svar models either to impose block
%   exogeneity or to impose other forms of linear restrictions. When not
%   empty, **estim_linear_restrictions** must be a 2-column cell:
%   - Each row of the first column represents a particular linear
%   combination of the estimated parameters. Those linear combinations are
%   constructed using the **coef** class. Check help for coef.coef for more
%   details.
%   - Each row of the second column holds the value of the linear
%   combination.
%
% - **estim_blocks** [{[]}|cell]: When not empty, this triggers blockwise
%   optimization. For further information on how to set blocks, see help
%   for dsge.create_estimation_blocks
%
% - **estim_priors** [{[]}|struct]: This provides an alternative to
%   setting priors inside the rise/dsge model file. Each field of the
%   structure must be the name of an estimated parameter. Each field will
%   hold a cell array whose structure is described in help
%   rise_generic.setup_priors.
%
% - **estim_penalty** [numeric|{1e+8}]: value of the objective function
%   when a problem occurs. Possible problems include:
%   - no solution found
%   - very low likelihood
%   - stochastic singularity
%   - problems computing the initial covariance matrix
%   - non-positive definite covariance matrices
%   - etc.
%
% - **estim_penalty_factor** [numeric|{10}]: when general nonlinear
% restrictions are present, RISE uses an estimation strategy in which the
% objective function is penalized as
% f_final=fval+estim_penalty_factor*sum(max(0,g)^2) where g is a vector of
% the values of the restrictions, which are expected to be of the form
% g(x)<=0. See **estim_general_restrictions** above.
%
% - **optimset** [struct]: identical to matlab's optimset
%
% - **optimizer** [char|function handle|cell]: This can be the name of a
%   standard matlab optimizer or RISE optimization routine or a
%   user-defined optimization procedure available of the matlab search
%   path. If the optimzer is provided as a cell, then the first element of
%   the cell is the name of the optimizer or its handle and the remaining
%   entries in the cell are additional input arguments to the user-defined
%   optimization routine. A user-defined optimization function should have
%   the following syntax ::
%      [xfinal,ffinal,exitflag]=optimizer(fh,x0,lb,ub,options,varargin);
%   That is, it accepts as inputs:
%       - **fh**: the function to optimize
%       - **x0**: a vector column of initial values of the parameters
%       - **lb**: a vector column of lower bounds
%       - **ub**: a vector column of upper bounds
%       - **options**: a structure of options whose fields will be similar
%           to matlab's optimset
%       - **varargin**: additional arguments to the user-defined
%           optimization procedure
%   That is, it provides as outputs:
%       - **xfinal**: the vector of final values
%       - **ffinal**: the value of **fh** at **xfinal**
%       - **exitflag**: a flag similar to the ones provided by matlab's
%       optimization functions.
%
% - **estim_optimizer_hessian** [{false}|true]: Use the hessian computed by
% the optimizer. Else, store the hessian returned by the optimizer, but
% also compute and use the numerical hessian.
%
% Outputs
% --------
%
% - **obj** [rise|dsge|rfvar|svar]: model object parameterized with the
%   mode found and holding additional estimation results and statistics
%   that can be found under obj.estimation
%
% More About
% ------------
%
% - recursive estimation may be done easily by passing a different
%   estim_end_date at the beginning of each estimation run.
%
% - It is also possible to estimate a dsge model using conditional
%   future information on endogenous (**forecast_cond_endo_vars**) as well.
%   as on exogenous (**forecast_cond_exo_vars**). The dataset provided in
%   this case must have several pages. The first page is the actual data,
%   while the subsequent pages are the expectations data.
%   See help dsge.forecast for more information
%
% Examples
% ---------
%
% See also: RISE_GENERIC/ESTIMATE

if ~isempty(obj)
    % Initially set the filtering/smoothing flag to false (during estimation).
    % This is especially important given that the objective function could be
    % optimal_simple_rule_posterior, in which case there is no filtering going
    % on.
    obj=set(obj,'kf_filtering_level',0);
end

obj=estimate@rise_generic(obj,varargin{:});

end
