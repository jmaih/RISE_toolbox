function obj=estimate(obj,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% varargin: data,optimset,optimizer,hessian_type,estim_start_from_mode,estim_parallel
% recursive estimation may be done easily by passing a different
% estim_end_date at the beginning of each estimation run.

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
        'estim_penalty',1e+8);
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

% make the hessian positive definite if necessary
if obj(1).options.hessian_repair
    [Hc,Hinv] = utils.hessian.conditioner(H);
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
    
    % capture the final parameters and their standard deviations
    obj(ii).estimation.posterior_maximization.mode=x1;
    obj(ii).estimation.posterior_maximization.mode_stdev=SD;
    
    % rebuild the parameter object. This can be done now coz estimation_under_way is set to false
    % log_mdd=.5*npar*log(2*pi)-.5*log(det(H))+log_post;
    % or alternatively
    log_mdd=.5*npar*log(2*pi)+.5*log(det(Hinv))+log_post;
    
    obj(ii).estimation.posterior_maximization.log_prior=log_prior(1);
    obj(ii).estimation.posterior_maximization.log_endog_prior=log_prior(2);
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


