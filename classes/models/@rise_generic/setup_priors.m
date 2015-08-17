function obj=setup_priors(obj,MyPriors,error_control)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also:

% different ways of setting priors outside the model file
% P.pname={start_value,lower_bound,upper_bound};
% P.pname={start_value,prior_mean,prior_stdev,'distribution'};
% P.pname={start_value,prior_mean,prior_stdev,'distribution',lower_bound};
% P.pname={start_value,prior_mean,prior_stdev,'distribution',lower_bound,upper_bound};
% P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)'};
% P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)',lower_bound};
% P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)',lower_bound,upper_bound};
%
% distribution can be any of the following: beta_pdf, cauchy_pdf,
% gamma_pdf, inv_gamma_pdf, laplace_pdf, left_triang_pdf, logistic_pdf,
% lognormal_pdf, normal_pdf, pareto_pdf, right_triang_pdf, uniform_pdf,
% weibull_pdf

if nargin<3
    error_control=[];
end
warnstate=warning('query','all');
warning('off','optim:fmincon:SwitchingToMediumScale')% %
warning('off','optimlib:fmincon:WillRunDiffAlg')
warning('off','optimlib:fmincon:SwitchingToMediumScaleBecauseNoGrad')

if ~isempty(fieldnames(MyPriors))
    
    disp(' ')
    disp('Now computing the hyperparameters for estimation...')
    disp(' ')
    
    param_tex_names=obj.parameters.tex_name;
    
    fields=fieldnames(MyPriors);
    priors=[];
    obj.estimation_restrictions=[];
    name_file_line=[];
    error_control_flag=~isempty(error_control);
    if error_control_flag
        error_control=[fields(:),error_control];
    end
    
    [estnames,is_dirichlet,dirichlet,error_control]=parameter_list(fields,...
        MyPriors,error_control);
    % linking estimated parameters to parameters
    %-------------------------------------------
    obj.estimation_restrictions=parameters_links(obj,estnames);
    new_dirichlet=utils.distrib.dirichlet_shortcuts();
    est_id=0;
    while est_id<numel(estnames)
        est_id=est_id+1;
        if error_control_flag
            name_file_line=error_control(est_id,:);
        end
        if is_dirichlet(est_id)
            % find the corresponding dirichlet and do all its elements and
            % increment est_id. Take the first dirichlet, use it and
            % destroy it.
            do_the_dirichlet()
        else
            do_one_typical()
        end
    end
    
    % for efficiency, this should be done at estimation time?...
    if ~isempty(priors)
        % load the distributions
        tmp={priors.prior_distrib};
        if ~isempty(tmp)
            effective_distributions=unique(tmp);
            distr_locs=cell(1,numel(effective_distributions));
            for ii=1:numel(effective_distributions)
                distr_locs{ii}=find(strcmp(effective_distributions{ii},tmp));
                % get the handle on the distributions but not for the
                % dirichlet: they need to be processed separately.
                if ~strcmp(effective_distributions{ii},'dirichlet')
                    lndens=distributions.(effective_distributions{ii})();
                    effective_distributions{ii}=lndens;
                end
            end
            obj.estim_hyperparams=[[priors.a]',[priors.b]'];
            obj.estim_distributions=effective_distributions;
            obj.estim_distrib_locations=distr_locs;
        end
        obj.estim_dirichlet=new_dirichlet;
    end
    
    warning(warnstate)
    
    obj.estimation=orderfields(...
        struct('endogenous_priors',[],'priors',{priors},...
        'posterior_maximization',struct(...
        'estim_start_time',[],'estim_end_time',[],...
        'log_lik',[],'log_post',[],'log_prior',[],'log_endog_prior',[],...
        'active_inequalities_number',0,'hessian',[],'vcov',[],'mode',[],...
        'mode_stdev',[],'funevals',[],...
        'log_marginal_data_density_laplace',[]...
        ),...
        'posterior_simulation',[])...
        );
end

    function do_the_dirichlet()
        d1=dirichlet(1);
        est_id=est_id-1;
        prior_trunc=obj.options.prior_trunc;
        for ii_=1:d1.n_1
            fildname=d1.names{ii_};
            [position,~,pname,chain,state]=decompose_parameter_name(obj,fildname,est_id==1);
            if isempty(position)
                error([fildname,' is not recognized as a parameter'])
            end
            [~,~,icdfn]=distributions.beta();
            bounds=[icdfn(prior_trunc,d1.a(ii_),d1.b(ii_)),...
                icdfn(1-prior_trunc,d1.a(ii_),d1.b(ii_))];
            est_id=est_id+1;
            priors(est_id)=struct('name',pname,'chain',chain,'state',state,...
                'start',d1.moments.mean(ii_),'lower_quantile',nan,...
                'upper_quantile',nan,'prior_mean',d1.moments.mean(ii_),...
                'prior_stdev',d1.moments.sd(ii_),...
                'prior_distrib','dirichlet',...
                'prior_prob',1,'lower_bound',bounds(1),...
                'upper_bound',bounds(2),...
                'tex_name',param_tex_names{position},'id',est_id,...
                'prior_trunc',prior_trunc,'a',d1.a(ii_),'b',d1.b(ii_));
            
            priors(est_id)=format_estimated_parameter_names(priors(est_id),param_tex_names{position});
            
            disp([' parameter: ',upper(pname),', density:',upper('dirichlet'),...
                ', hyperparameters: [',num2str(d1.a(ii_)),' ',num2str(d1.b(ii_)),'],',...
                'convergence ',num2str(0)])
        end
        new_dirichlet(end+1)=utils.distrib.dirichlet_shortcuts(dirichlet(1).a,...
            dirichlet(1).location);
        dirichlet=dirichlet(2:end);
    end

    function do_one_typical()
        fildname=estnames{est_id};
        tmp=MyPriors.(fildname);
        n_entries=numel(tmp);
        if ~iscell(tmp)||numel(tmp)<3
            error('all fields of the prior structure should be cell arrays with at least 3 elements')
        end
        [position,~,pname,chain,state]=decompose_parameter_name(obj,fildname,est_id==1);
        if isempty(position)
            error([fildname,' is not recognized as a parameter'])
        end
        start=parser.push_if_validated(tmp{1},@(x)isfinite(x),'start value',name_file_line);
        lq=nan;    uq=nan;     pmean=nan;    pstdev=nan;    prior_prob=1;
        distrib='uniform';    lb=nan;    ub=nan;
        if n_entries==3
            lb=parser.push_if_validated(tmp{2},@(x)isfinite(x),'lower bound',name_file_line);
            ub=parser.push_if_validated(tmp{3},@(x)isfinite(x)&& x>lb,'upper bound value',name_file_line);
            lq=lb;
            uq=ub;
        else
            distrib=parser.push_if_validated(tmp{4},@(x)ischar(x),'distribution(prob)',name_file_line);
            left_par=strfind(distrib,'(');
            if isempty(left_par)
                left_par=length(distrib)+1;
                pmean=parser.push_if_validated(tmp{2},@(x)isfinite(x),'prior mean',name_file_line);
                pstdev=parser.push_if_validated(tmp{3},@(x)isfinite(x)&& x>0,'prior stdev',name_file_line);
                % we have to change default for the probability
                prior_prob=nan;
            else
                right_par=strfind(distrib,')');
                prior_prob=eval(distrib(left_par+1:right_par-1));
                prior_prob=parser.push_if_validated(prior_prob,@(x)isfinite(x) && x>=0 && x<=1,'prior probability',name_file_line);
                lq=parser.push_if_validated(tmp{2},@(x)isfinite(x),'lower quantile',name_file_line);
                uq=parser.push_if_validated(tmp{3},@(x)isfinite(x)&& x>lq,'upper quantile',name_file_line);
            end
            distrib=distrib(1:left_par-1);
            if n_entries>4
                lb=parser.push_if_validated(tmp{5},@(x)isfinite(x),'lower bound',name_file_line);
                if n_entries>5
                    ub=parser.push_if_validated(tmp{6},@(x)isfinite(x),'upper bound',name_file_line);
                    if n_entries>6
                        error('number of entries in setting up the prior cannot exceed 6')
                    end
                end
            end
        end
        block=struct('name',pname,'chain',chain,'state',state,'start',start,...
            'lower_quantile',lq,'upper_quantile',uq,'prior_mean',pmean,...
            'prior_stdev',pstdev,'prior_distrib',distrib,'prior_prob',prior_prob,...
            'lower_bound',lb,'upper_bound',ub);
        % Adapt the name before setting the prior
        %----------------------------------------
        block=format_estimated_parameter_names(block,param_tex_names{position});
        % hyperparameters and other things
        %--------------------------------
        priors=prior_setting_engine(priors,block,est_id,obj.options.prior_trunc);
    end
end

function block=format_estimated_parameter_names(block,par_tex_name)
block.tex_name=par_tex_name;
if ~strcmp(block.chain,'const')
    RegimeState_tex=['(',block.chain,',',sprintf('%0.0f',block.state),')'];
    RegimeState=['_',block.chain,'_',sprintf('%0.0f',block.state)];
    block.name=[block.name,RegimeState];
    block.tex_name=strcat(block.tex_name,' ',RegimeState_tex);
end
end

function [estnames,is_dirichlet,dirichlet,error_control]=parameter_list(...
    fields,MyPriors,error_control)

is_dirichlet=strncmp(fields,'dirichlet',9);
fields=fields(:).';
dirichlet=struct('a',{},'b',{},'moments',{},'n_1',{},...
    'pointers',{},'location',{},'names',{});
ndirich=sum(is_dirichlet);
if ndirich
    error_control_flag=~isempty(error_control);
    name_count=0;
    n=1000;
    estnames=cell(1,n);
    if error_control_flag
        tmp=error_control;
        ncols=size(tmp,2);
        error_control=cell(n,ncols);
    end
    dirich_count=0;
    for icol=1:numel(fields)
        dname=fields{icol};
        if is_dirichlet(icol)
            dirich_count=dirich_count+1;
            vals=MyPriors.(dname);
            s02=vals{1}.^2;
            vals=reshape(vals(2:end),2,[]);
            pnames=vals(1,:);
            m_main=cell2mat(vals(2,:));
            m0=1-sum(m_main);
            a_sum=m0*(1-m0)/s02-1;
            if a_sum<=0
                error(['dirichlet # ',int2str(dirich_count),...
                    ' appears to have a too big standard deviation'])
            end
            m=[m_main(:);m0];
            a=a_sum*m;
            [dirichlet(dirich_count).a,dirichlet(dirich_count).b,...
                dirichlet(dirich_count).moments,...
                dirichlet(dirich_count).fval,...
                dirichlet(dirich_count).space]=distributions.dirichlet(a);
            dirichlet(dirich_count).pointers=1:numel(pnames);
            dirichlet(dirich_count).n_1=numel(m)-1;
            dirichlet(dirich_count).names=pnames;
        else
            pnames={dname};
        end
        n_names=numel(pnames);
        if name_count+n_names>=n
            estnames{n+100}={};
            n=n+100;
        end
        pos=name_count+(1:n_names);
        if error_control_flag
            item=tmp(icol*ones(1,n_names),:);
            if is_dirichlet(icol)
                item(:,1)=pnames(:);
            end
            error_control(pos,:)=item;
        end
        estnames(pos)=pnames;
        name_count=pos(end);
        if is_dirichlet(icol)
            dirichlet(dirich_count).location=pos;
        end
    end
    estnames=estnames(1:name_count);
    error_control=error_control(1:name_count,:);
    is_dirichlet=false(1,name_count);
    is_dirichlet([dirichlet.location])=true;
else
    estnames=fields;
end
end

function prior=prior_setting_engine(prior,parray,id,prior_trunc)
if nargin<5
    prior_trunc=1e-10;
end

% for truncation
invgamma_upper_bound_truncation=10;

parray.id=id;
parray.prior_distrib=strrep(parray.prior_distrib,'_pdf','');
parray.prior_trunc=prior_trunc;
mean_stdev_flag=isnan(parray.prior_prob);
if mean_stdev_flag
    lqtl_mean=parray.prior_mean;
    uqtl_std=parray.prior_stdev;
else
    if parray.prior_prob<=0||parray.prior_prob>1
        error([mfilename,':: probability for parameter ',parray.name,' should be in (0,1]'])
    end
    lqtl_mean=parray.lower_quantile;
    uqtl_std=parray.upper_quantile;
end

% find the hyperparameters
[parray.a,parray.b,moments,ffinal]=distributions.(parray.prior_distrib)(lqtl_mean,uqtl_std,parray.prior_prob);
parray.prior_mean=moments.mean;
parray.prior_stdev=moments.sd;

disp([' parameter: ',upper(parray.name),', density:',upper(parray.prior_distrib),...
    ', hyperparameters: [',num2str(parray.a),' ',num2str(parray.b),'],',...
    'convergence ',num2str(ffinal)])
% get the functions of the distribution
[~,~,icdfn]=distributions.(parray.prior_distrib)();
bounds=[icdfn(prior_trunc,parray.a,parray.b),icdfn(1-prior_trunc,parray.a,parray.b)];
if isempty(parray.prior_mean)||isempty(parray.prior_stdev)
    disp([mfilename,'(GENTLE WARNING):: for these hyperparameters, the distribution ',...
        'does not have well-defined moments'])
end
the_message='';
if ismember(parray.prior_distrib,{'inv_gamma'})
    if bounds(2)>invgamma_upper_bound_truncation
        the_message=[mfilename,'(GENTLE WARNING):: upper bound of inverse gamma distribution ',...
            'truncated at ',num2str(invgamma_upper_bound_truncation)];
    end
    bounds(2) = min(bounds(2),invgamma_upper_bound_truncation);
end
if isfinite(parray.lower_bound),bounds(1)=max(bounds(1),parray.lower_bound);end
if isfinite(parray.upper_bound),bounds(2)=min(bounds(2),parray.upper_bound);end
% if the distribution has been truncated, say it here.
disp(the_message)

% check that the starting value is not outside the bounds
%--------------------------------------------------------
if any(parray.start<bounds(1))||any(parray.start>bounds(2))
    error([mfilename,':: parameter ',parray.name,' (',num2str(parray.start),') outside its bounds [',num2str(bounds),']'])
end
parray.lower_bound = bounds(1);
parray.upper_bound = bounds(2);

if isempty(prior)
    prior=parray;
else
    prior(end+1)=parray;
end

end
