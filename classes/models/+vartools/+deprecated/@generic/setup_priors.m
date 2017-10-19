function obj=setup_priors(obj,MyPriors,error_control)
% SETUP_PRIORS -- format the priors for the estimation of models in RISE.
%
% Syntax
% -------
% ::
%
%   obj=SETUP_PRIORS(obj,MyPriors)
%
%   obj=SETUP_PRIORS(obj,MyPriors,error_control)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object. This function is not
% meant to be called directly by the user. But it shows how to write the
% priors when those are not written inside the model file.
%
% - **MyPriors** [struct]: Each field of the structure is one of the
% following alternatives
%   - Maximum likelihood or uniform priors
%       P.pname={start_value,lower_bound,upper_bound};
%   - Bayesian prior using mean and standard deviation
%       P.pname={start_value,prior_mean,prior_stdev,'distribution'};
%   - Same as above except for adding a hard lower bound
%       P.pname={start_value,prior_mean,prior_stdev,'distribution',lower_bound};
%   - Same as above except for adding a hard upper bound
%       P.pname={start_value,prior_mean,prior_stdev,'distribution',lower_bound,upper_bound};
%   - Bayesian prior using quantiles of the distribution
%       P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)'};
%   - Same as above except for adding a hard lower bound
%       P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)',lower_bound};
%   - Same as above except for adding a hard upper bound
%       P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)',lower_bound,upper_bound};
%
% - **error_control** [empty|cell]: element constructed by RISE for
% controling the syntax used in the model file.
%
% Outputs
% --------
%
% - **obj** [rise|dsge|svar|rfvar]: updated model object.
%
% More About
% ------------
%
% - This function is also indirectly used for svar and rfvar objects.
%
% Examples
% ---------
%
% - Possible distributions include: beta, cauchy, gamma, inv_gamma, laplace,
% left_triang, logistic, lognormal, normal, pareto, right_triang,
% uniform, weibull
%
% See also: RISE_GENERIC/ESTIMATE, DSGE/ESTIMATE

if nargin<3
    
    error_control=[];
    
end

nobj=numel(obj);

if nobj>1
    
    for iobj=1:nobj
        
        obj(iobj)=setup_priors(obj(iobj),MyPriors,error_control);
        
    end
    
    return
    
end

warnstate=warning('query','all');
warning('off','optim:fmincon:SwitchingToMediumScale')% %
warning('off','optimlib:fmincon:WillRunDiffAlg')
warning('off','optimlib:fmincon:SwitchingToMediumScaleBecauseNoGrad')

fields=fieldnames(MyPriors);

if isempty(fields)
    
    return
    
end

disp(' ')
disp('Now computing the hyperparameters for estimation...')
disp(' ')

param_tex_names=obj.parameters.tex_name;

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

while est_id < numel(estnames)
    
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

obj.estim_priors_data=utils.prior.load_priors(struct(),priors,new_dirichlet);

warning(warnstate)

obj.estimation.priors=priors;

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
            % We set 0 and 1 as lower and upper bounds
            %------------------------------------------
            bounds=[prior_trunc,1-prior_trunc];
            est_id=est_id+1;
            mean_i=d1.moments.mean(ii_);
            est_struct=struct('name',pname,'chain',chain,'state',state,...
                'start',mean_i,'lower_quantile',nan,...
                'upper_quantile',nan,'prior_mean',mean_i,...
                'prior_stdev',d1.moments.sd(ii_),...
                'prior_distrib','dirichlet',...
                'prior_prob',1,'lower_bound',bounds(1),...
                'upper_bound',bounds(2),...
                'tex_name',param_tex_names{position},'id',est_id,...
                'prior_trunc',prior_trunc,'a',d1.a(ii_),'b',d1.b(ii_));
            if est_id==1
                priors=est_struct;
            else
                priors(est_id)=est_struct;
            end
            priors(est_id)=format_estimated_parameter_names(priors(est_id),param_tex_names{position});
            
            disp([' parameter: ',upper(pname),', density:',upper('dirichlet'),...
                ', hyperparameters: [',num2str(d1.a(ii_)),' ',num2str(d1.b(ii_)),'],',...
                'convergence ',num2str(0)])
        end
        new_dirichlet(end+1)=utils.distrib.dirichlet_shortcuts(dirichlet(1).a,...
            dirichlet(1).location,[],[],dirichlet(1).s);
        dirichlet=dirichlet(2:end);
    end

    function do_one_typical()
        
        fildname=estnames{est_id};
        
        tmp=MyPriors.(fildname);
        
        [position,~,pname,chain,state]=decompose_parameter_name(obj,fildname,est_id==1);
        
        if isempty(position)
            
            error([fildname,' is not recognized as a parameter'])
            
        end
        
        block=utils.prior.cell2block(tmp,pname,chain,state,name_file_line);
        % Adapt the name before setting the prior
        %----------------------------------------
        block=format_estimated_parameter_names(block,param_tex_names{position});
        % hyperparameters and other things
        %--------------------------------
        priors=utils.prior.prior_setting_engine(priors,block,est_id,obj.options.prior_trunc);
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
    'pointers',{},'location',{},'names',{},'s',{});

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
            
            h=numel(m);
            
            a=a_sum*m;
            
            [dirichlet(dirich_count).a,dirichlet(dirich_count).b,...
                dirichlet(dirich_count).moments,...
                dirichlet(dirich_count).fval,...
                dirichlet(dirich_count).space]=distributions.dirichlet(a);
            
            dirichlet(dirich_count).pointers=1:numel(pnames);
            
            dirichlet(dirich_count).n_1=h-1;
            
            dirichlet(dirich_count).names=pnames;
            
            % Un-normalized weight of the diagonal term.
            dirichlet(dirich_count).s=utils.distrib.dirichlet_sum_weights(h);
            
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
    
    if error_control_flag
        
        error_control=error_control(1:name_count,:);
        
    end
    
    is_dirichlet=false(1,name_count);
    
    is_dirichlet([dirichlet.location])=true;
    
else
    
    estnames=fields;
    
end

end

