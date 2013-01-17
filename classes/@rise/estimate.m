function obj=estimate(obj,varargin)
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
        'estim_start_time',[],...
        'estim_end_time',[],...
        'estim_mode_file',[],...
        'estim_general_restrictions',[],... % holds a function that takes as input the model object and returns the
        'estim_blocks',[]);
    % violations of the restrictions in a vector for instance in order to impose the max operator ZLB, Linde and Maih
    obj=mergestructures(estimation_engine(),...
        main_defaults);
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
                obj(jj).options.optimset=mysetfield(obj(jj).options.optimset,varargin{2*ii});
            end
            discard=[2*ii-1,2*ii];
        end
    end
    varargin(discard)=[];
    if ~isempty(varargin)
        for jj=1:nobj
            obj(jj).options=mysetfield(obj(jj).options,varargin{:});
        end
    end
end

estim_start_time=clock;
for jj=1:nobj
    if isempty(obj(jj).options.estim_start_time)
        obj(jj).options.estim_start_time=estim_start_time;
    end
end

optim_options=obj(1).options.optimset;
optimizer=obj(1).options.optimizer;
hessian_type=obj(1).options.hessian_type;
estim_start_from_mode=obj(1).options.estim_start_from_mode;
estim_parallel=obj(1).options.estim_parallel;
estim_general_restrictions=obj(1).options.estim_general_restrictions;
estim_blocks=obj(1).options.estim_blocks;
if ~isempty(estim_blocks)
    estim_blocks=create_estimation_blocks(obj(1),estim_blocks);
end

% Load the function that computes the likelihood
% preliminary checks:
param_names={obj(1).estimated_parameters.name};
% Load the function that computes the likelihood
for ii=1:nobj
    if ~isequal(param_names,{obj(ii).estimated_parameters.name})
        error([mfilename,':: optimization parameters should be the same across all models and ordered in the same way'])
    end
    % this will be useful when deciding to check for stability in the markov switching model
    obj(ii).estimation_under_way=true;
    % %     isunif=strcmp('uniform',{obj(ii).estimated_parameters.distribution});
    % %     is_uniform=[is_uniform,isunif(:)];
    % load the mode
    obj(ii)=load_mode(obj(ii));
    % Initially set the filtering/smoothing flag to false (during estimation).
    % This is especially important given that the objective function could be
    % optimal_simple_rule_posterior, in which case there is no filtering going
    % on.
    obj(ii).options.kf_filtering_level=0;
    if ~isempty(obj(ii).options.estim_general_restrictions)
        % collect the information about the degree of filtering
        % 0 (no filters), 1(filtered), 2(filtered+updated), 3(filtered+updated+smoothed)
        obj(ii).options.kf_filtering_level=obj(ii).options.estim_general_restrictions();
        if obj(ii).options.kf_filtering_level && obj(ii).is_optimal_simple_rule_model
            error([mfilename,':: Cannot do filtering under estimation of optimal simple rules'])
        end
    end
    % Now destroy the parameter object
    obj(ii).parameters=object2cell(obj(ii).parameters);
end

% this will record the different problems encounter during estimation
list_of_issues=cell(0);

if ~isempty(estim_start_from_mode) && ~islogical(estim_start_from_mode)
    estim_start_from_mode=logical(estim_start_from_mode);
end
% prior_trunc=obj.options.prior_trunc;
% debug=obj(1).options.debug;

[obj,issue,retcode]=load_data(obj);
if retcode
    if ~all([obj.is_optimal_simple_rule_model])
        error([mfilename,':: ',decipher_error(retcode)])
    end
end
if ~isempty(issue)
    list_of_issues=[list_of_issues;{issue}];
end
% eventually I should put options for recursive estimation
% and for passing new data to the model on the fly.

xmode=vertcat(obj(1).estimated_parameters.mode);
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
        x0=vertcat(obj(1).estimated_parameters.startval);
    otherwise
end

lb=vertcat(obj(1).estimated_parameters.lb);
ub=vertcat(obj(1).estimated_parameters.ub);

% This does not work well under parallel mode, especially for old versions
% of Matlab
try %#ok<TRYNC>
    prior_plots(obj)
end

[nonlcon,nconst]=reprocess_nolinear_restrictions(obj(1).parameter_random_draws_restrictions);
% the functions in obj(1).parameter_random_draws_restrictions only check
% whether the restrictions are violated or not but do not give the strength
% of the violation, which we need in order to apply DEB. So we need to use
% the reprocessed form.

npar=size(x0,1);
Nsim=max(1,estim_parallel);
x0=[x0,nan(npar,Nsim-1)];
f0=nan(1,Nsim);
% all objects are evaluated at the same point. Without a second argument,
% this is exactly what will happen.
[~,f0(1),~,retcode0,viol]=big_wrapper(x0(:,1));
if retcode0||any(viol)
    % first check constraint violations
    f0(1)=obj(1).options.Penalty;
end
% the objects are automatically updated and potentially contain crucial
% information going forward. In particular, they contain information
% about whether the models are stationary or not.

if f0(1)<obj(1).options.Penalty
    beg=2;
else
    beg=1;
end

%% disable those elements
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:illConditionedMatrix')

fprintf(1,'%s\n','Looking for good enough start values. Please wait...');
for ii=beg:Nsim
    NotDone=true;
    iter=0;
    while NotDone
        [xtest,ftest,retcode]=generate_starting_point(@big_wrapper);
        if ftest<obj(1).options.Penalty
            NotDone=false;
            f0(ii)=ftest;
            x0(:,ii)=xtest;
        end
        iter=iter+1;
        if iter>=obj(1).options.estim_max_trials
            error([mfilename,':: No admissible starting value after ',...
                int2str(obj(1).options.estim_max_trials),' trials'])
        else
            fprintf(1,'%3.0d :: %s\n',iter,decipher_error(retcode));
        end
    end
    disp(['Starting value # ',int2str(ii),' found after ',int2str(iter),' iterations'])
    ratio=ii/Nsim;
    fprintf(1,'%s\n',['...', num2str(100*ratio),'% done']);
end


[x1,f1,H,issue,viol,obj]=big_wrapper(x0,'estimate'); %#ok<ASGLU>
viol=viol(viol>0);
numberOfActiveInequalities=numel(viol);

% make the hessian positive definite if necessary
[Hc,Hinv] = hessian_conditioner(H);
if ~isequal(Hc,H)
    warning([mfilename,':: non-positive definite hessian made diagonally dominant']) %#ok<WNTAG>
    H=Hc; clear Hc
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
    obj(ii).options.estim_end_time=estim_end_time;
    % set the filtering/smoothing flag to 3 in order to get out the final
    % outputs with the smoothed and filtered variables (if feasible)
    if ~obj(ii).is_optimal_simple_rule_model
        obj(ii).options.kf_filtering_level=3;
    end
    [log_post,log_lik,log_prior,~,~,obj(ii)]=log_posterior_kernel(obj(ii),x1);
    
    if obj(ii).options.kf_filtering_level
        % put the filters in the time series format
        obj(ii)=save_filters(obj(ii));
    end
    % capture the final parameters and their standard deviations
    for jj=1:npar
        obj(ii).estimated_parameters(jj)=set_properties(obj(ii).estimated_parameters(jj),...
            'mode',x1(jj),'mode_std',SD(jj));
    end
    % rebuild the parameter object
    obj(ii).parameters=cell2object(obj(ii).parameters,'rise_param');
    % log_mdd=.5*npar*log(2*pi)-.5*log(det(H))+log_post;
    % or alternatively
    log_mdd=.5*npar*log(2*pi)+.5*log(det(Hinv))+log_post;
        
    obj(ii)=obj(ii).set_properties('log_prior',log_prior(1),...
        'log_endog_prior',log_prior(2),...
        'log_post',log_post,...
        'log_lik',log_lik,...
        'log_mdd_laplace',log_mdd,...
        'numberOfActiveInequalities',numberOfActiveInequalities,...
        'Hessian',H,...
        'vcov',Hinv);
    
    obj(ii).list_of_issues=list_of_issues;
    
    % now we change the flag so that stability can be tested
    obj(ii).estimation_under_way=false;
    
    save([obj(ii).options.results_folder,filesep,'estimation',filesep,...
        'estimated_model'],'obj','x1','x0','f1','f0','H','Hinv')
    
end
% disp Estimation results
print_estimation_results(obj);

warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:illConditionedMatrix')

    function [x1,f1,H,issue,viol,objLast]=big_wrapper(x0,action)
        if nargin<2
            action='eval';
            if nargin<1
                x0=lb+(ub-lb).*rand(npar,1);
            end
        end
        violLast=[];
        objLast=[];
        retcodeLast=[];
        viol=[];
        switch action
            case 'estimate'
                %         [x,f,eflag,output]=problem.optimizer(@fh_wrapper,@nonlcon_with_gradient);
                PROBLEM_=struct('objective',@fh_wrapper,...
                    'x0',x0,...
                    'lb',lb,...
                    'ub',ub,...
                    'nonlcon',@nonlcon_with_gradient,... % the nonlinear constraints restrictions take the same inputs as fh_wrapper
                    'options',optim_options,...
                    'solver',optimizer);
                
                [x1,f1,H,issue]=estimation_engine(PROBLEM_,hessian_type,estim_blocks);
            case 'eval'
                [f1,retcode_0]=fh_wrapper(x0);
                x1=x0;
                H=[];
                issue=retcode_0;
                viol=nonlcon_with_gradient(x0);
            otherwise
                error(['unknown type of action ',action])
        end
        
        function [minus_log_post,retcode]=fh_wrapper(x)
            % this function returns the minimax if there are many objects
            % evaluated simultaneously
            
            % all objects are evaluated at the same point
            % the reason you want to output the object here is because it potentially
            % contains crucial information going forward. In particular, it contains
            % information about whether the model is stationary or not.
            fval=obj(1).options.Penalty*ones(1,nobj);
            for mo=1:nobj
                [fval(mo),~,~,~,retcode,obj(mo)]=log_posterior_kernel(obj(mo),x);
                if retcode
                    break
                end
            end
            % Now take the negative for minimization
            minus_log_post=-min(fval);
            objLast=obj;
            retcodeLast=retcode;
            % nonlinear constraints might be incompatible with blockwise
            % optimization. In that case, it is better to compute the
            % restriction violation while evaluating the objective and save
            % the results to pass on to the optimizer when it calls the
            % nonlinear constraints.
            violLast=nonlcon_with_gradient(x,true);
        end
        
        function [viol,grad]=nonlcon_with_gradient(x,do_it)
            grad=[];
            if nargin<2
                do_it=false;
            end
            if do_it
                if retcodeLast
                    viol=ones(nconst,1)*realmax/nconst;
                    % this should make it hard enough for this vector to stay in any population
                else
                    viol=nonlcon(x);
                end
                if ~isempty(estim_general_restrictions)
                    viols=estim_general_restrictions(objLast);
                    viol=[viol(:);viols(:)];
                    % it may be important to know the number of
                    % restrictions returned by this function or at least
                    % the user should return the proper number of arguments
                    % whenever something fails.
                end
            else
                viol=violLast;
            end
        end
    end
end

    function [nonlcon,nconst]=reprocess_nolinear_restrictions(nonlcon)
        % transform the nonlinear constraints. I would like to keep the
        % flexibility of knowning what parameters enter the constraints and
        % so I do not do this in format parameters
        nonlcon=nonlcon(:,2);
        nconst=size(nonlcon,1);
        for iconstr=1:nconst
            nonlcon{iconstr}=strrep(func2str(nonlcon{iconstr}),'@(param)','');
            % now remove inequalities
            cutoff_type={'>=','<=','>','<','='};
            for itype=1:numel(cutoff_type)
                cutoff_locs=strfind(nonlcon{iconstr},cutoff_type{itype});
                if ~isempty(cutoff_locs)
                    cutoff_type=cutoff_type{itype};
                    break
                end
            end
            if ~isempty(cutoff_locs)
                span=length(cutoff_type);
                left=nonlcon{iconstr}(1:cutoff_locs-1);
                right=nonlcon{iconstr}(cutoff_locs+span:end);
                switch cutoff_type
                    case '>='
                        nonlcon{iconstr}=[right,'-(',left,')-eps;'];
                    case '<='
                        nonlcon{iconstr}=[left,'-(',right,')-eps;'];
                    case '>'
                        nonlcon{iconstr}=[right,'-(',left,');'];
                    case '<'
                        nonlcon{iconstr}=[left,'-(',right,');'];
                    case '='
                        nonlcon{iconstr}=['abs(',left,'-(',right,'))-eps;'];
                end
            end
        end
        nonlcon=cell2mat(nonlcon(:)');
        nonlcon=nonlcon(1:end-1);
        if isempty(nonlcon)
            nonlcon='0';
            nconst=1;
        end
        nonlcon=str2func(['@(param)[',nonlcon,']']);
    end
% function rst=push_time_series(C)
% rst=rise_time_series(C{1},C{2},C{3});
% end

