function obj=setup_priors(obj,MyPriors,error_control)% is_switching=[ParameterInfo.is_switching];
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
    
    param_names=obj.parameters.name;
    param_tex_names=obj.parameters.tex_name;
    regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));
    chain_names=obj.markov_chains.chain_names;
    par_nbr=sum(obj.parameters.number);
    
    fields=fieldnames(MyPriors);
    priors=[];
    obj.estimation_restrictions=[];
    name_file_line=[];
    error_control_flag=~isempty(error_control);
    if error_control_flag
        error_control=[fields(:),error_control];
    end
    for est_id=1:numel(fields)
        tmp=MyPriors.(fields{est_id});
        if error_control_flag
            name_file_line=error_control(est_id,:);
        end
        n_entries=numel(tmp);
        if ~iscell(tmp)||numel(tmp)<3
            error('all fields of the prior structure should be cell arrays with at least 3 elements')
        end
        [position,regime_states,pname,chain,state]=decompose_name(fields{est_id});
        if isempty(position)
            error([fields{est_id},' is not recognized as a parameter'])
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
        % linking estimated parameters to parameters
        %-------------------------------------------
        obj.estimation_restrictions=[obj.estimation_restrictions
            position+(regime_states(:)-1)*par_nbr,est_id*ones(numel(regime_states),1)];
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
                % get the handle on the distributions
                lndens=distributions.(effective_distributions{ii})();
                effective_distributions{ii}=lndens;
            end
            obj.estim_hyperparams=[[priors.a]',[priors.b]'];
            obj.estim_distributions=effective_distributions;
            obj.estim_distrib_locations=distr_locs;
        end
    end
    
    warning(warnstate)
    
    obj.estimation=orderfields(struct('endogenous_priors',[],...
        'estim_start_time',[],'estim_end_time',[],...
        'log_lik',[],'log_post',[],'log_prior',[],'log_endog_prior',[],...
        'log_marginal_data_density',struct('laplace',[],'modified_harmonic_mean',[],'chib_jeliazkov',[]),...
        'active_inequalities_number',0,'priors',{priors},...
        'hessian',[],'vcov',[],'mode',[],'mode_stdev',[],...
        'funevals',[]));
end

    function block=format_estimated_parameter_names(block,par_tex_name)
        block.tex_name=par_tex_name;
        if ~strcmp(block.chain,'const')
            RegimeState=['(',block.chain,',',sprintf('%0.0f',block.state),')'];
            block.name=[block.name,RegimeState];
            block.tex_name=strcat(block.tex_name,' ',RegimeState);
        end
    end

    function [position,regime_states,pname,chain,state]=decompose_name(pname)
        position=find(strcmp(pname,param_names));
        if ~isempty(position)
            state=1;
            chain='const';
            chain_id=find(strcmp(chain,chain_names));
        else
            ptex=parser.valid_param_name_to_tex_name(pname,chain_names);
            left_par=strfind(ptex,'(');
            right_par=strfind(ptex,')');
            pname=ptex(1:left_par-1);
            comma=strfind(ptex,',');
            position=find(strcmp(pname,param_names));
            if isempty(left_par)||isempty(right_par)||isempty(comma)||isempty(position)
                error(['"',ptex(1:left_par-1),'" is not recognized as a parameter name'])
            end
            chain=ptex(left_par+1:comma-1);
            state=str2double(ptex(comma+1:right_par-1));
            chain_id=find(strcmp(chain,chain_names));
        end
        govChain=obj.parameters.governing_chain(position);
        if ~(chain_id==govChain)
            error(['parameter ',pname,' is not controlled by ',chain_names(chain_id)])
        end
        % locate the state in the regimes
        regime_states=find(regimes(:,chain_id)==state);
        if isempty(regime_states)
            error([sprintf('%0.0f',state),' is not a valid state for parameter ',pname])
        end
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
% if mean_stdev_flag
%     parray.prior_prob=0.9;
%     bounds=distributions.find_bounds(parray.prior_distrib,...
%         parray.prior_mean,...
%         parray.prior_stdev,...
%         parray.prior_prob);
%     parray.lower_quantile=bounds(1);
%     parray.upper_quantile=bounds(2);
% end
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
    error([mfilename,':: parameter ',name,' (',num2str(parray.start),') outside its bounds [',num2str(bounds),']'])
end
parray.lower_bound = bounds(1);
parray.upper_bound = bounds(2);

if isempty(prior)
    prior=parray;
else
    prior(end+1)=parray;
end

end
