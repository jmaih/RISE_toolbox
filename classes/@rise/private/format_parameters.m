function obj=format_parameters(obj,ParameterInfo,ParameterizationArray,...
    MarkovChains,Param_rest_block)

% re-order the names of the markov chains right when creating it...

param_names={ParameterInfo.name};
param_tex_names={ParameterInfo.tex_name};
is_switching=[ParameterInfo.is_switching];
obj.markov_chains=struct('name',MarkovChains(1,:));
n_states=cell2mat(MarkovChains(2,:));

[obj.Regimes,obj.journal]=chain_grid(n_states);
[reg_nbr,chain_nbr]=size(obj.Regimes);

for ch=1:chain_nbr
    % each markov chain has states
    for jj=1:n_states(ch)
        % and each state of the markov chain may appear in several regimes
        obj.markov_chains(ch).states(jj).regime=find(obj.Regimes(:,ch)==jj);
    end
end

par_nbr=numel(param_names);
obj.parameters=rise_param.empty(0,0);
obj.estimated_parameters=rise_estim_param.empty(0,0);

model_is_parameterized=~isempty(ParameterizationArray);
if model_is_parameterized
    ParameterizationArray=cell2struct(ParameterizationArray,...
        {'names','chains','states','startval','lb','ub','distr','probs'},2);
    estimated=~isnan([ParameterizationArray.lb])&~isnan([ParameterizationArray.ub]);
    cumsum_estimated=cumsum(estimated);
end

warnstate=warning('query','all');                
warning('off','optim:fmincon:SwitchingToMediumScale')% %
warning('off','optimlib:fmincon:WillRunDiffAlg')
warning('off','optimlib:fmincon:SwitchingToMediumScaleBecauseNoGrad')

obj.estimation_restrictions=[];
estim_iter=0;
NotParameterized={};
for ii=1:par_nbr
    pname=param_names{ii};
    p_tex_name=param_tex_names{ii};
    value=nan(1,reg_nbr);
    if model_is_parameterized
        locs=find(strcmp(pname,{ParameterizationArray.names}));
        if isempty(locs)
            NotParameterized=[NotParameterized,{pname}]; %#ok<AGROW>
        else
            chain_name=ParameterizationArray(locs(1)).chains;
            val_i=[ParameterizationArray(locs).startval];
            plb_i=[ParameterizationArray(locs).lb];
            pub_i=[ParameterizationArray(locs).ub];
            prob_i=[ParameterizationArray(locs).probs];
            states_i=[ParameterizationArray(locs).states];
            chain_id=locate_variables(chain_name,{obj.markov_chains.name});
            for jj=1:n_states(chain_id)
                % the chain is in state jj
                chain_state=jj;
                % in which regimes is the state in chain_state?
                lhs=obj.markov_chains(chain_id).states(chain_state).regime;
                % now locate the position of state chain_state in the vector of
                % states
                rhs=states_i==chain_state;
                % now put the particular values of that states in the different
                % regimes on the left-hand side
                value(lhs)=val_i(rhs);
                if ~isnan(plb_i(rhs)) && ~isnan(pub_i(rhs)) % then the parameter is estimated
                    estim_iter=estim_iter+1;
                    est_distr=ParameterizationArray(locs(rhs)).distr;
                    [formated_par_name,formated_par_tex_name]=...
                        format_estimated_parameter_names(pname,p_tex_name,chain_name,chain_state);
                    
                    est_id=cumsum_estimated(locs(rhs));
                    obj.estimated_parameters(estim_iter,1)=...
                        rise_estim_param(formated_par_name,formated_par_tex_name,est_id,...
                        val_i(rhs),plb_i(rhs),pub_i(rhs),est_distr,prob_i(rhs),obj.options.prior_trunc);
                    obj.estimation_restrictions=[obj.estimation_restrictions
                        ii+(lhs-1)*par_nbr,est_id*ones(numel(lhs),1)];
                end
            end
        end
    end
    obj.parameters(ii,1)=rise_param('name',pname,'tex_name',p_tex_name,'id',ii,'startval',value,'is_switching',is_switching(ii));
end
warning(warnstate) 

% identification restrictions
n_restr=numel(Param_rest_block);
direct=cell(n_restr,2);
estim_param_names={obj.estimated_parameters.name};
delete_equation_in_draws=false(1,numel(Param_rest_block));
for ii=1:n_restr
    equation_i=Param_rest_block{ii};
    equation_d=Param_rest_block{ii};
    for jj=1:size(equation_i,2)
        if isempty(equation_i{2,jj})
            p_id=find(strcmp(equation_i{1,jj},estim_param_names));
            if ~isempty(p_id)
                p_id=obj.estimated_parameters(p_id).id;
                equation_d{1,jj}=['param(',int2str(p_id),')'];
                % gathering the indexes of the parameters involved
                direct{ii,1}=[direct{ii,1},p_id];
            end  
        else
            pname=equation_i{1,jj};
            state=equation_i{2,jj}{2};
            chain_=equation_i{2,jj}{1};
            p_id=find(strcmp(pname,param_names));
            c_id= strcmp(chain_,MarkovChains(1,:));
            col=find(obj.Regimes(:,c_id)==state,1,'first');
            % Now we just need to replace the correct parameter location in
            % the matrix...
            equation_i{1,jj}=['M(',int2str(p_id),',',int2str(col),')'];
            % location in the estimated parameters
            if strcmp(chain_,'const')
                p_id=find(strcmp(pname,estim_param_names));
            else
                p_id=find(strcmp([pname,'(',chain_,',',int2str(state),')'],estim_param_names));
            end
            if isempty(p_id)
                delete_equation_in_draws(ii)=true;
            else
                p_id=obj.estimated_parameters(p_id).id;
                equation_d{1,jj}=['param(',int2str(p_id),')'];
                % gathering the indexes of the parameters involved
                direct{ii,1}=[direct{ii,1},p_id];
            end            
        end
    end
    direct{ii,1}=unique(direct{ii,1});
    Param_rest_block{ii}=cell2mat(equation_i(1,:));
    direct{ii,2}=cell2mat(equation_d(1,:));
    direct{ii,2}=str2func(['@(param)',direct{ii,2}]);
end
% throw out the restrictions that involve non estimated parameters as we
% won't be able to check that they hold or not when drawing random vectors
% of parameters.
direct(delete_equation_in_draws,:)=[];
obj.parameter_restrictions=Param_rest_block;
obj.parameter_random_draws_restrictions=direct;
% check that all the parameters in the model are in use
not_in_use=param_names(~obj.is_in_use_parameter);
if ~isempty(not_in_use)
    warning([mfilename,' :: the following parameters do not seem to affect model ',obj.filename,...
        '. You may want to discard them from model for tidiness'])
    disp(not_in_use)
end
if ~isempty(NotParameterized)
    warning([mfilename,' :: the following parameters are not assigned values in model ',obj.filename,...
        '. The model might not solve'])  %#ok<*WNTAG>
    disp(NotParameterized)
end
in_use_but_not_parameterized=intersect(param_names(obj.is_in_use_parameter),NotParameterized);
if ~isempty(in_use_but_not_parameterized)
    disp(in_use_but_not_parameterized)
    warning([mfilename,upper([':: The parameters above affect model ',obj.filename,' but are not assigned values. I am very angry!!!'])])
end

% now resort the estimated parameter vector using their ids
obj.estimated_parameters([obj.estimated_parameters.id])=obj.estimated_parameters;

% load the distributions
tmp={obj.estimated_parameters.distribution};
if ~isempty(tmp)
    effective_distributions=unique(tmp);
    distr_locs=cell(1,numel(effective_distributions));
    for ii=1:numel(effective_distributions)
        distr_locs{ii}=find(strcmp(effective_distributions{ii},tmp));
        % get the handle on the distributions
        lndens=distributions.(effective_distributions{ii})();
        effective_distributions{ii}=lndens;
    end
    
    hyperparams=nan(numel(obj.estimated_parameters),4);
    hyperparams(:,1)=[obj.estimated_parameters.a];
    hyperparams(:,2)=[obj.estimated_parameters.b];
    hyperparams(:,3)=[obj.estimated_parameters.c];
    hyperparams(:,4)=[obj.estimated_parameters.d];
    obj.estim_hyperparams=hyperparams;
    obj.estim_distributions=effective_distributions;
    obj.estim_distrib_locations=distr_locs;
end

% measurement errors restrictions: this has to come after the
% parameters
obj.measurement_errors_restrictions=[];
for ii=1:obj.NumberOfObservables(1) % pick only the endogenous observables
    vi=obj.varobs(ii).name;
    loc=find(strcmp(['stderr_',vi],{obj.parameters.name}));
    % the line above is a bit inefficient. I have to change it
    if ~isempty(loc)
        obj.measurement_errors_restrictions=...
            [obj.measurement_errors_restrictions;ii,loc];
    end
end


function [formated_par_name,formated_par_tex_name]=...
    format_estimated_parameter_names(par_name,par_tex_name,chain_name,chain_state)

% now construct the rise_estim_param objects
formated_par_name=par_name;
formated_par_tex_name=par_tex_name;
if ~strcmp(chain_name,'const')
    RegimeState=['(',chain_name,',',int2str(chain_state),')'];
    formated_par_name=[formated_par_name,RegimeState];
    formated_par_tex_name=strcat(formated_par_tex_name,' ',RegimeState);
end

