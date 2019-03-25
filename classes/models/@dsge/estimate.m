function [obj,filtration]=estimate(obj,varargin)
% Estimate the parameters of a RISE model
%
% ::
%
%   obj=estimate(obj)
%   obj=estimate(obj,varargin)
%   [obj,filtration]=estimate(...)
%
% Args:
%
%    obj (rise | dsge | rfvar | svar): model object
%    varargin : additional optional inputs among which the most relevant
%      for estimation are:
%
%       - **estim_parallel** [integer | {1}]: Number of starting values
%
%       - **estim_start_from_mode** [true | false | {[]}]: when empty, the user is
%         prompted to answer the question as to whether to start estimation from
%         a previously found mode or not. If true or false, no question is asked.
%
%       - **estim_start_date** [numeric | char | serial date]: date of the first
%         observation to use in the dataset provided for estimation
%
%       - **estim_end_date** [numeric | char | serial date]: date of the last
%         observation to use in the dataset provided for estimation
%
%       - **estim_max_trials** [integer | {500}]: When the initial value of the
%         log-likelihood is too low, RISE uniformly draws from the prior support
%         in search for a better starting point. It will try this for a maximum
%         number of **estim_max_trials** times before squeaking with an error.
%
%       - **estim_start_vals** [{[]} | struct]: when not empty, the parameters
%         whose names are fields of the structure will see their start values
%         updated or overriden by the information in **estim_start_vals**. There
%         is no need to provide values to update the start values for the
%         estimated parameters.
%
%       - **estim_general_restrictions** [{[]} | function handle | cell array]: when
%         not empty, the argument can be a function handle or a cell array
%         containing the function handle and additional input arguments. The
%         general syntax for the calling the function handle is
%         viol=myfunc(obj,varargin), with **obj** the parameterized RISE object
%         which will be used in the computation of the restrictions violations.
%         Hence the restrictions are entered either as  @myfunc or as
%         {@myfunc,arg2,arg3,...}. Originally, RISE will call the function
%         without any inputs. In that case, RISE expects the output to be a
%         structure with fields :
%
%          - **number_of_restrictions** : number of restrictions
%          - **kf_filtering_level** [0 | 1 | 2 | 3]: if 0, no filters are required
%            for the computation of the restrictions. If 1, only the filtered
%            variables are required. If 2, the updated variables are required.
%            If 3, the smoothed variables are required.
%
%         When the function is called with inputs, RISE expects as output the
%         values of the restrictions. The sign of the violations does not matter.
%         All the user has to do is to put a zero where the restrictions are not
%         violated.
%
%       - **estim_linear_restrictions** [{[]} | cell]: This is most often used in
%         the estimation of rfvar or svar models either to impose block
%         exogeneity or to impose other forms of linear restrictions. When not
%         empty, **estim_linear_restrictions** must be a 2-column cell:
%
%          - Each row of the first column represents a particular linear
%            combination of the estimated parameters.
%          - Each row of the second column holds the value of the linear
%            combination.
%
%       - **estim_nonlinear_restrictions** [{[]} | cell]: When not
%         empty, **estim_nonlinear_restrictions** must be a k x 1 cell, with each
%         row representing a particular restriction on the parameters. e.g. for a
%         switching model, one can have alpha(zlb,1)>alpha(zlb,2), which can also
%         be written as alpha_zlb_1>alpha_zlb_2.
%         The restrictions can also be equality restrictions. In this case,
%         however, it is assumed that the parameters entering the lhs of
%         restrictions are not estimated. e.g. alpha(zlb,1)=3*cos(alpha(zlb,2))+1.
%
%       - **estim_endogenous_priors** [{[]} | function handle]: When not empty,
%         **estim_endogenous_priors** must be a function handle such that when
%         called without inputs, it returns a struct with fields:
%
%       - **estim_priors** [{[]}|struct]: This provides an alternative to
%         setting priors inside the rise/dsge model file. Each field of the
%         structure must be the name of an estimated parameter. Each field will
%         hold a cell array whose structure is described in help
%         GENERIC/setup_priors.
%
%          - **priors** : cell array of estimation priors. more explicitly, each
%            entry of the cell array is itself a cell with the same syntax as the
%            priors for estimation, EXCEPT the start value!
%          - **kf_filtering_level** [0 | 1 | 2 | 3]: if 0, no filters are required
%            for the computation of the endogenous priors. If 1, only the
%            filtered variables are required. If 2, the updated variables are
%            required. If 3, the smoothed variables are required.
%
%         When the function handle is called with TWO inputs, what is
%         returned is a vector of values for which RISE will evaluate the
%         endogenous  prior. This vector should have the same length as the
%         previous cell array. The two inputs are:
%         - **obj** : The model object
%         - **filtration** : a structure containing the filters
%
%       - **estim_blocks** [{[]} | cell]: When not empty, this triggers blockwise
%         optimization. For further information on how to set blocks, see help
%         for dsge.create_estimation_blocks
%
%       - **estim_penalty** [numeric | {1e+8}]: value of the objective function
%         when a problem occurs. Possible problems include:
%
%          - no solution found
%          - very low likelihood
%          - stochastic singularity
%          - problems computing the initial covariance matrix
%          - non-positive definite covariance matrices
%          - etc.
%
%       - **estim_penalty_factor** [numeric | {10}]: when general nonlinear
%         restrictions are present, RISE uses an estimation strategy in which the
%         objective function is penalized as
%         f_final=fval+estim_penalty_factor*sum(max(0,g)^2) where g is a vector
%         of the values of the restrictions, which are expected to be of the form
%         g(x)<=0. See **estim_general_restrictions** above.
%
%       - **optimset** [struct]: identical to matlab's optimset
%
%       - **optimizer** [char | function handle | cell|{fmincon}]: This can
%         be the name of a standard matlab optimizer or RISE optimization
%         routine or a user-defined optimization procedure available of the
%         matlab search path. If the optimzer is provided as a cell, then
%         the first element of the cell is the name of the optimizer or its
%         handle and the remaining entries in the cell are additional input
%         arguments to the user-defined optimization routine. A
%         user-defined optimization function should have the following
%         syntax :: 
%
%            [xfinal,ffinal,exitflag,H]=optimizer(fh,x0,lb,ub,options,varargin);
%
%         That is, it accepts as inputs:
%
%             - **fh**: the function to optimize
%             - **x0**: a vector column of initial values of the parameters
%             - **lb**: a vector column of lower bounds
%             - **ub**: a vector column of upper bounds
%             - **options**: a structure of options whose fields will be similar
%               to matlab's optimset
%             - **varargin**: additional arguments to the user-defined
%               optimization procedure
%
%         That is, it provides as outputs:
%
%             - **xfinal**: the vector of final values
%             - **ffinal**: the value of **fh** at **xfinal**
%             - **exitflag**: a flag similar to the ones provided by matlab's
%               optimization functions.
%             - **H**: an estimate of the Hessian
%
%       - **estim_barrier** [{false} | true]: never allow constraints to be
%         violated in no circumstances.
%
% Returns:
%    :
%
%    - **obj** [rise | dsge | rfvar | svar]: model object parameterized with the
%      mode found and holding additional estimation results and statistics
%      that can be found under obj.estimation
%
%    - **filtration** [struct]: structure with the filtration information
%      of the model parameterized at the mode
%
% Note:
%
%    - recursive estimation may be done easily by passing a different
%      estim_end_date at the beginning of each estimation run.
%
%    - It is also possible to estimate a dsge model using conditional
%      future information on endogenous (**forecast_cond_endo_vars**) as well.
%      as on exogenous (**forecast_cond_exo_vars**). The dataset provided in
%      this case must have several pages. The first page is the actual data,
%      while the subsequent pages are the expectations data.
%      See help dsge.forecast for more information
%
%
% See also: GENERIC/SETUP_PRIORS

if isempty(obj)
    
    mydefaults=estimate@generic(obj,varargin{:});
    
    mydefaults=[mydefaults
        {'estim_priors',[],@(x)isstruct(x),...
        'estim_priors must be a structure'}];
        
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj
        
        disp_defaults(mydefaults);
        
    end

    
    return
    
else
    % Initially set the filtering/smoothing flag to false (during estimation).
    % This is especially important given that the objective function could be
    % optimal_simple_rule_posterior, in which case there is no filtering going
    % on.
    obj=set(obj,'kf_filtering_level',0);
    
    [obj,filtration]=estimate@generic(obj,varargin{:});
    
end


end
