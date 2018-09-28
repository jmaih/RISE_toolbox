function [obj,filtration]=estimate(obj,varargin)

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
    
    filtration=cell(1,nobj);
    
    if NumWorkers
        
        pctRunOnAll('utils.estim.warnings_disable();')
        
        pctRunOnAll('utils.estim.prior.warnings_fmincon_disable();')
        
    end
    
    parfor (ii=1:nobj,NumWorkers)
        
        [obj(ii),filtration{ii}]=estimate(obj(ii),varargin{:}); %#ok<PFBNS>
        
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

filtration=cell(1,nobj);

for ii=1:nobj
    % set the filtering/smoothing flag to 3 in order to get out the final
    % outputs with the smoothed and filtered variables (if feasible)
    
    [log_post,log_lik,log_prior,~,~,obj(ii),xmode,filtration{ii}]=conclude_estimation(obj(ii),x1);
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
    
    post_max.x0=x0;
    
    post_max.f0=f0;
    
    obj(ii).estimation.posterior_maximization=post_max;
    
    obj(ii).list_of_issues=list_of_issues;
    
end

if nobj==1
    
    filtration=filtration{1};
    
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
    }; %#ok<ISCLSTR>

end
