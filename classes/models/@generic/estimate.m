function obj=estimate(obj,varargin)
% ESTIMATE - estimates the parameters of a RISE model
%
% ::
%
%
%   obj=ESTIMATE(obj)
%   obj=ESTIMATE(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: model object
%
%    - **varargin** additional optional inputs among which the most relevant
%      for estimation are:
%
%    - **estim_parallel** [integer|{1}]: Number of starting values
%
%    - **estim_start_from_mode** [true|false|{[]}]: when empty, the user is
%      prompted to answer the question as to whether to start estimation from
%      a previously found mode or not. If true or false, no question is asked.
%
%    - **estim_start_date** [numeric|char|serial date]: date of the first
%      observation to use in the dataset provided for estimation
%
%    - **estim_end_date** [numeric|char|serial date]: date of the last
%      observation to use in the dataset provided for estimation
%
%    - **estim_max_trials** [integer|{500}]: When the initial value of the
%      log-likelihood is too low, RISE uniformly draws from the prior support
%      in search for a better starting point. It will try this for a maximum
%      number of **estim_max_trials** times before squeaking with an error.
%
%    - **estim_start_vals** [{[]}|struct]: when not empty, the parameters
%      whose names are fields of the structure will see their start values
%      updated or overriden by the information in **estim_start_vals**. There
%      is no need to provide values to update the start values for the
%      estimated parameters.
%
%    - **estim_general_restrictions** [{[]}|function handle|cell array]: when
%      not empty, the argument can be a function handle or a cell array
%      containing the function handle and additional input arguments. The
%      general syntax for the calling the function handle is
%      viol=myfunc(obj,varargin), with **obj** the parameterized RISE object
%      which will be used in the computation of the restrictions violations.
%      Hence the restrictions are entered either as  @myfunc or as
%      {@myfunc,arg2,arg3,...}. Originally, RISE will call the function
%      without any inputs. In that case, RISE expects the output to be a
%      structure with fields :
%          - **number_of_restrictions** : number of restrictions
%          - **kf_filtering_level** [0|1|2|3]: if 0, no filters are required
%          for the computation of the restrictions. If 1, only the filtered
%          variables are required. If 2, the updated variables are required.
%          If 3, the smoothed variables are required.
%      When the function is called with inputs, RISE expects as output the
%      values of the restrictions. The sign of the violations does not matter.
%      All the user has to do is to put a zero where the restrictions are not
%      violated.
%
%    - **estim_linear_restrictions** [{[]}|cell]: This is most often used in
%      the estimation of rfvar or svar models either to impose block
%      exogeneity or to impose other forms of linear restrictions. When not
%      empty, **estim_linear_restrictions** must be a 2-column cell:
%      - Each row of the first column represents a particular linear
%      combination of the estimated parameters.
%      - Each row of the second column holds the value of the linear
%      combination.
%
%    - **estim_nonlinear_restrictions** [{[]}|cell]: When not
%      empty, **estim_nonlinear_restrictions** must be a k x 1 cell, with each
%      row representing a particular restriction on the parameters. e.g. for a
%      switching model, one can have alpha(zlb,1)>alpha(zlb,2), which can also
%      be written as alpha_zlb_1>alpha_zlb_2.
%      The restrictions can also be equality restrictions. In this case,
%      however, it is assumed that the parameters entering the lhs of
%      restrictions are not estimated. e.g. alpha(zlb,1)=3*cos(alpha(zlb,2))+1.
%
%    - **estim_endogenous_priors** [{[]}|function handle]: When not empty,
%      **estim_endogenous_priors** must be a function handle such that when
%      called without inputs, it returns a struct with fields:
%      - **priors** : cell array of estimation priors. more explicitly, each
%      entry of the cell array is itself a cell with the same syntax as the
%      priors for estimation, EXCEPT the start value!
%      - **kf_filtering_level** [0|1|2|3]: if 0, no filters are required
%          for the computation of the endogenous priors. If 1, only the
%          filtered variables are required. If 2, the updated variables are
%          required. If 3, the smoothed variables are required.
%      When the function handle is called with an input, what is returned is a
%      vector of values for which RISE will evaluate the endogenous prior.
%      This vector should have the same length as the previous cell array.
%
%    - **estim_blocks** [{[]}|cell]: When not empty, this triggers blockwise
%      optimization. For further information on how to set blocks, see help
%      for dsge.create_estimation_blocks
%
%    - **estim_penalty** [numeric|{1e+8}]: value of the objective function
%      when a problem occurs. Possible problems include:
%      - no solution found
%      - very low likelihood
%      - stochastic singularity
%      - problems computing the initial covariance matrix
%      - non-positive definite covariance matrices
%      - etc.
%
%    - **estim_penalty_factor** [numeric|{10}]: when general nonlinear
%      restrictions are present, RISE uses an estimation strategy in which the
%      objective function is penalized as
%      f_final=fval+estim_penalty_factor*sum(max(0,g)^2) where g is a vector
%      of the values of the restrictions, which are expected to be of the form
%      g(x)<=0. See **estim_general_restrictions** above.
%
%    - **optimset** [struct]: identical to matlab's optimset
%
%    - **optimizer** [char|function handle|cell]: This can be the name of a
%      standard matlab optimizer or RISE optimization routine or a
%      user-defined optimization procedure available of the matlab search
%      path. If the optimzer is provided as a cell, then the first element of
%      the cell is the name of the optimizer or its handle and the remaining
%      entries in the cell are additional input arguments to the user-defined
%      optimization routine. A user-defined optimization function should have
%      the following syntax ::
%         [xfinal,ffinal,exitflag]=optimizer(fh,x0,lb,ub,options,varargin);
%      That is, it accepts as inputs:
%          - **fh**: the function to optimize
%          - **x0**: a vector column of initial values of the parameters
%          - **lb**: a vector column of lower bounds
%          - **ub**: a vector column of upper bounds
%          - **options**: a structure of options whose fields will be similar
%              to matlab's optimset
%          - **varargin**: additional arguments to the user-defined
%              optimization procedure
%      That is, it provides as outputs:
%          - **xfinal**: the vector of final values
%          - **ffinal**: the value of **fh** at **xfinal**
%          - **exitflag**: a flag similar to the ones provided by matlab's
%          optimization functions.
%
%    - **estim_barrier** [{false}|true]: never allow constraints to be
%    violated in no circumstances.
%
% Returns:
%    :
%
%    - **obj** [rise|dsge|rfvar|svar]: model object parameterized with the
%      mode found and holding additional estimation results and statistics
%      that can be found under obj.estimation
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
% Example:
%
%    See also:

nobj=numel(obj);

if nobj==0 % same as isempty(obj)
    
    mydefaults=the_defaults();
    
    mydefaults=[mydefaults
        optimization.estimation_engine()];
    %         load_data(obj),...
    
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
elseif nobj>1
    
    NumWorkers=utils.parallel.get_number_of_workers();
    
    parfor (ii=1:nobj,NumWorkers)
        
        obj(ii)=estimate(obj(ii),varargin{:}); %#ok<PFBNS>
        
    end
    
    return
    
end

nn=length(varargin);

if nn
    
    if rem(nn,2)~=0
        
        error([mfilename,':: arguments must come in pairs'])
        
    end
    
    if ~isempty(varargin)
        
        for jj=1:nobj
            
            obj(jj)=set(obj(jj),varargin{:});
            
        end
        
    end
    
end

estim_start_time=clock;

for jj=1:nobj
    
    if isempty(obj(jj).estimation.posterior_maximization.estim_start_time)
        % maybe we used several optimizers one after the other. In that
        % case the start of the estimation should not change
        obj(jj).estimation.posterior_maximization.estim_start_time=estim_start_time;
        
    end
    
end

estim_start_from_mode=obj(1).options.estim_start_from_mode;

% Load the function that computes the likelihood
% preliminary checks:
param_names={obj(1).estimation.priors.name};
% Load the function that computes the likelihood
for ii=1:nobj
    
    if ~isequal(param_names,{obj(ii).estimation.priors.name})
        
        error([mfilename,':: optimization parameters should be the same across all models and ordered in the same way'])
    
    end
    % this will be useful when deciding to check for stability in the markov switching model
    obj(ii).estimation_under_way=true;
    % load the mode
    obj(ii)=load_mode(obj(ii));

end

% this will record the different problems encounter during estimation
list_of_issues=cell(0);

if ~isempty(estim_start_from_mode) && ~islogical(estim_start_from_mode)
    
    estim_start_from_mode=logical(estim_start_from_mode);
    
end

[obj,issue,retcode]=load_data(obj);

if retcode
    
    if ~all([obj.is_optimal_simple_rule_model])
        
        error([mfilename,':: ',utils.error.decipher(retcode)])
    
    end
    
end

if ~isempty(issue)
    
    list_of_issues=[list_of_issues;{issue}];
    
end

xmode=obj(1).estimation.posterior_maximization.mode;

response='n';

if ~isempty(xmode)
    
    if ~isempty(estim_start_from_mode)
        
        if estim_start_from_mode
            
            response='y';
            
        else
            
            response='n';
            
        end
        
    else
        
        response='';
        
        while ~ismember(response,{'y','n'})
            
            response=input('previous estimation found. Do you want to load the mode? [y/n]','s');
        
        end
        
    end
    
end

switch response
    
    case 'y'
        
        x0=xmode;
    
    case 'n'
        
        x0=[obj(1).estimation.priors.start];
    
    otherwise
        
end

x0=x0(:);

lb=[obj(1).estimation.priors.lower_bound]; lb=lb(:);

ub=[obj(1).estimation.priors.upper_bound]; ub=ub(:);

% transform initial conditions
%-----------------------------
[obj,x0,lb_short,ub_short]=transform_parameters(obj,x0,lb,ub);

% do posterior maximization
%--------------------------
[x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb_short,ub_short); %#ok<ASGLU>

% extend the output but not the short hessian
%---------------------------------------------
linear_restricts=obj(1).linear_restrictions_data;

% x1 = linear_restricts.a_func(x1);
% 
% x0 = linear_restricts.a_func(x0); %#ok<NASGU>
% H=linear_restricts.a_func(H,true);

x1=unstransform_parameters(obj(1),x1);

x0=unstransform_parameters(obj(1),x0); %#ok<NASGU>

numberOfActiveInequalities=numel(viol);

if ~isempty(issue)
    
    list_of_issues=[list_of_issues;{issue}];

end

estim_end_time=clock;

for ii=1:nobj
    % set the filtering/smoothing flag to 3 in order to get out the final
    % outputs with the smoothed and filtered variables (if feasible)
    
    [log_post,log_lik,log_prior,~,~,obj(ii),xmode]=conclude_estimation(obj(ii),x1);
    % compute the penalties for the restrictions violations
    %-------------------------------------------------------
    g=evaluate_general_restrictions(obj(ii));
    
    nonlin_penalty=utils.estim.penalize_violations(g{1},obj(ii).options.estim_penalty_factor);
    
    post_max=obj(ii).estimation.posterior_maximization;
    % set the end time for estimation, even if it does not include the
    % smoothing part done in conclude_estimation.
    post_max.estim_end_time=estim_end_time;
    % capture the final parameters and their standard deviations
    post_max.mode=xmode;
    
    % Hessian
    post_max.hessian=H;
    
    post_max.log_prior=log_prior(1);
    post_max.log_endog_prior=log_prior(2);
    post_max.nonlinear_restrictions_penalty=nonlin_penalty;
    post_max.log_post=log_post;
    post_max.log_lik=log_lik;
    post_max.active_inequalities_number=numberOfActiveInequalities;
    post_max.funevals=funevals;
    
    % Marginal data density and other variable outputs that can be calculated separately
    %-----------------------------------------------------------------------------------
    post_max=generic_tools.posterior_maximization_variable_quantities(...
        post_max,linear_restricts.a_func);
    
    obj(ii).estimation.posterior_maximization=post_max;
    
    obj(ii).list_of_issues=list_of_issues;
    
    if isdir(obj(ii).options.results_folder)
        
        save([obj(ii).options.results_folder,filesep,'estimation',filesep,...
            'estimated_model'],'obj','x1','x0','f1','f0','H')
    
    end
    
end
% disp Estimation results
print_estimation_results(obj);

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'estim_barrier',false,@(x)islogical(x),...
    'estim_barrier must be a logical'
    
    'estim_start_from_mode',[],@(x)islogical(x),...
    'estim_start_from_mode must be a logical'
    
    'estim_penalty_factor',10,@(x)num_fin(x) && x >0,...
    'estim_penalty_factor must be a finite and positive scalar'
    
    'estim_penalty',1e+8,@(x)num_fin(x) && x >0,...
    'estim_penalty must be a finite and positive scalar'
    
    'estim_max_trials',500,@(x)num_fin_int(x) && x >0,...
    'estim_max_trials must be a finite and positive integer'
    
    'estim_parallel',1,@(x)num_fin_int(x) && x>0,...
    'estim_parallel must be a positive integer' % number of starting values
    
    'estim_start_date','',@(x)isempty(x)||is_date(x)||is_serial(x),...
    'estim_start_date must be a valid date'
    
    'estim_end_date','',@(x)isempty(x)||is_date(x)||is_serial(x),...
    'estim_end_date must be a valid date'
    
    'estim_start_vals',[],@(x)isstruct(x),...
    'estim_start_vals must be a struct'
    
    'estim_general_restrictions',[],...
    @(x)ischar(x)||iscell(x)||isa(x,'function_handle'),...
    'estim_general_restrictions must be a char, a cell or a function handle' % holds a function that takes as input the model object and returns the
    
    'estim_linear_restrictions',[],@(x)iscell(x)&&size(x,2)==2,...
    'estim_linear_restrictions must be a two-column cell'
    
    'estim_nonlinear_restrictions',[],@(x)iscellstr(x),...
    'estim_nonlinear_restrictions must be a cell array of strings'
    
    'estim_blocks',[],@(x)iscell(x),...
    'estim_blocks must be a cell array with groupings of parameter names'
    
    'estim_endogenous_priors',[],@(x)isa(x,'function_handle'),...
    'estim_endogenous_priors must be a functio handle'
    };

end
