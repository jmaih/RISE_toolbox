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
% - **hessian_type** [{'fd'}|'opg']: The hessian is either computed by
%   finite differences (fd) or by outer-product-gradient (opg)
%
% - **hessian_repair** [{false}|true]: If the Hessian is not positive
%   definite, it nevertheless can be repaired and prepared for a potential
%   mcmc simulation.
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
% See also: 


nobj=numel(obj);
if nobj==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    main_defaults=struct('estim_parallel',1,... % number of starting values
        'estim_start_from_mode',[],...
        'estim_start_date','',...
        'estim_end_date','',...
        'estim_max_trials',500,... % maximum number of trials when generating admissible starting values
        'estim_start_vals',[],...
        'estim_general_restrictions',[],... % holds a function that takes as input the model object and returns the
        'estim_linear_restrictions',[],...
        'estim_blocks',[],...
        'estim_priors',[],...
        'estim_penalty',1e+8,...
        'estim_penalty_factor',10);
    % violations of the restrictions in a vector for instance in order to impose the max operator ZLB, Linde and Maih
    obj=utils.miscellaneous.mergestructures(optimization.estimation_engine(),...
        main_defaults);
    %         load_data(obj),...
    return
end
nn=length(varargin);
if nn
    if rem(nn,2)~=0
        error([mfilename,':: arguments must come in pairs'])
    end
    discard=[];
    for ii=1:.5*nn
        if strcmp(varargin{2*ii-1},'optimset')
            for jj=1:nobj
                obj(jj).options.optimset=utils.miscellaneous.setfield(obj(jj).options.optimset,varargin{2*ii});
            end
            discard=[2*ii-1,2*ii];
        end
    end
    varargin(discard)=[];
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
    % %     isunif=strcmp('uniform',{obj(ii).estimation.priors.prior_distrib});
    % %     is_uniform=[is_uniform,isunif(:)];
    % load the mode
    obj(ii)=load_mode(obj(ii));
    % Initially set the filtering/smoothing flag to false (during estimation).
    % This is especially important given that the objective function could be
    % optimal_simple_rule_posterior, in which case there is no filtering going
    % on.
    obj(ii).options.kf_filtering_level=0;
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

obj=setup_linear_restrictions(obj);
obj=setup_general_restrictions(obj);

npar=size(x0,1);

[x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb,ub); %#ok<ASGLU>

numberOfActiveInequalities=numel(viol);

% make the hessian positive definite if necessary but keep a copy of the
% original
orig_H=H;
if obj(1).options.hessian_repair
    [Hc,Hinv] = utils.cov.conditioner(H);
    if max(abs(Hc(:)-H(:)))>1e-6
        warning([mfilename,':: non-positive definite hessian made diagonally dominant']) %#ok<WNTAG>
        H=Hc; clear Hc
        Hinv=inv(H);
    end
else
    Hinv=inv(H);
end
% the standard deviations
SD=sqrt(diag(Hinv));

if ~isempty(issue)
    list_of_issues=[list_of_issues;{issue}];
end

estim_end_time=clock;
for ii=1:nobj
    % set the end time for estimation, even if it does not include the
    % smoothing part below.
    obj(ii).estimation.posterior_maximization.estim_end_time=estim_end_time;
    % set the filtering/smoothing flag to 3 in order to get out the final
    % outputs with the smoothed and filtered variables (if feasible)
    
    [log_post,log_lik,log_prior,~,~,obj(ii)]=conclude_estimation(obj(ii),x1);
    % compute the penalties for the restrictions violations
    %-------------------------------------------------------
    g=evaluate_nonlinear_restrictions(obj(ii));
    nonlin_penalty=utils.estim.penalize_violations(g{1},obj(ii).options.estim_penalty_factor);
    
    % capture the final parameters and their standard deviations
    obj(ii).estimation.posterior_maximization.mode=x1;
    obj(ii).estimation.posterior_maximization.mode_stdev=SD;
    
    % rebuild the parameter object. This can be done now coz estimation_under_way is set to false
    % log_mdd=.5*npar*log(2*pi)-.5*log(det(H))+log_post;
    % or alternatively
    log_mdd=.5*npar*log(2*pi)+.5*log(det(Hinv))+log_post;
    
    obj(ii).estimation.posterior_maximization.log_prior=log_prior(1);
    obj(ii).estimation.posterior_maximization.log_endog_prior=log_prior(2);
    obj(ii).estimation.posterior_maximization.nonlinear_restrictions_penalty=nonlin_penalty;
    obj(ii).estimation.posterior_maximization.log_post=log_post;
    obj(ii).estimation.posterior_maximization.log_lik=log_lik;
    obj(ii).estimation.posterior_maximization.log_marginal_data_density_laplace=log_mdd;
    obj(ii).estimation.posterior_maximization.active_inequalities_number=numberOfActiveInequalities;
    obj(ii).estimation.posterior_maximization.vcov=Hinv;
    obj(ii).estimation.posterior_maximization.funevals=funevals;
    
    obj(ii).list_of_issues=list_of_issues;
    
    if isdir(obj(ii).options.results_folder)
        save([obj(ii).options.results_folder,filesep,'estimation',filesep,...
            'estimated_model'],'obj','x1','x0','f1','f0','H','Hinv')
    end
end
% disp Estimation results
print_estimation_results(obj);

end


