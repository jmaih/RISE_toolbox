function obj=load_functions(obj)

if ~exist(obj.options.results_folder,'dir')
    mkdir(obj.options.results_folder)
end
if ~exist([obj.options.results_folder,filesep,'graphs'],'dir')
    mkdir([obj.options.results_folder,filesep,'graphs'])
end
if ~exist([obj.options.results_folder,filesep,'estimation'],'dir')
    mkdir([obj.options.results_folder,filesep,'estimation'])
end
if ~exist([obj.options.results_folder,filesep,'simulations'],'dir')
    mkdir([obj.options.results_folder,filesep,'simulations'])
end

% Get rid of definitions
defcell=cell(0,2);
for ii=1:numel(obj.definitions)
    def_=obj.definitions(ii).shadow_dynamic;
    equality=strfind(def_,'=');
    % get rid of the semicolon
    defcell=[defcell;{def_(1:equality-1),def_(equality+1:end-1)}]; %#ok<AGROW>
end
for ii=1:numel(obj.definitions)
    for jj=ii+1:numel(obj.definitions)
        defcell{jj,2}=strrep(defcell{jj,2},defcell{ii,1},['(',defcell{ii,2},')']);
    end
end

% initialize
handle_struct=struct();

% function handles
if obj.is_dsge_var_model
    handle_struct.likelihood=@likelihood_dsge_var;
    obj.dsge_prior_weight_id=find(strcmp('dsge_prior_weight',{obj.parameters.name}),1);
elseif obj.is_optimal_simple_rule_model
    handle_struct.likelihood=@likelihood_optimal_simple_rule;
else
    handle_struct.likelihood=@likelihood_markov_switching_dsge;
end

%% parameter restrictions
% for the moment,I only allow the matrix of parameters. the idea is that if
% the restrictions are violated, evaluation should fail and the steady
% state should not be computed. In general, one could think of allowing the
% steady state to enter this game, but I already have enough problems with
% defining the steady state in markov switching, etc.
obj.parameter_restrictions=transpose(obj.parameter_restrictions(:));
tmp=['@(M)4*any([',cell2mat(transpose(obj.parameter_restrictions(:))),']==0)'];
% restrictions are violated if the function returns 4, which is the retcode
if ~isempty(obj.parameter_restrictions)
    tmp(isspace(tmp))=[];
    semcols=strfind(tmp,';');
    tmp(semcols(end))=[]; % remove last semicolon
    tmp=substitute_definitions(tmp,defcell);
end
handle_struct.parameter_restrictions=str2func(tmp);

%% static
tmp=['@(y,x,ss,param)[',cell2mat({obj.equations.shadow_static}),']'];
tmp(isspace(tmp))=[];
tmp(end-1)=[]; % remove last semicolon
tmp=substitute_definitions(tmp,defcell);
handle_struct.static=str2func(tmp);
%% balanced growth path
tmp=cellstr(obj.equations(1).shadow_balanced_growth_path)';
for ii=2:obj.NumberOfEquations
    tmp=[tmp,cellstr(obj.equations(ii).shadow_balanced_growth_path)']; %#ok<AGROW>
end
tmp=['@(y,x,ss,param)[',cell2mat(tmp),']'];
tmp(isspace(tmp))=[];
tmp(end-1)=[]; % remove last semicolon
tmp=substitute_definitions(tmp,defcell);
handle_struct.balanced_growth=str2func(tmp);

%% dynamic
tmp=['@(z,ss,param)[',cell2mat({obj.equations.shadow_dynamic}),']'];
tmp(isspace(tmp))=[];
tmp(end-1)=[]; % remove last semicolon
tmp=substitute_definitions(tmp,defcell);
% replace all x's
iter=nnz(obj.Lead_lag_incidence);
for ii=1:obj.NumberOfExogenous
    iter=iter+1;
    tmp=strrep(tmp,['x(',int2str(ii),')'],['y(',int2str(iter),')']);
end
% nvz=iter;
% replace all switching parameters
switching=find([obj.parameters.is_switching]);
for ii=1:numel(switching)
    iter=iter+1;
    tmp=strrep(tmp,['param(',int2str(switching(ii)),')'],['y(',int2str(iter),')']);
end
% replace y by z
tmp=strrep(tmp,'y(','z(');
% non vectorized form
handle_struct.dynamic=str2func(tmp);
% final steps for the vectorizations
for ii=1:iter
    tmp=strrep(tmp,['z(',int2str(ii),')'],['z(',int2str(ii),',:)']);
end
tmp=strrep(tmp,'^','.^');
tmp=strrep(tmp,'/','./');
tmp=strrep(tmp,'*','.*');
handle_struct.vectorized_dynamic=str2func(tmp);
%% transition matrix
handle_struct.transition_matrix=rise_anonymous.empty(0,1);
for ii=1:numel(obj.shadow_transition_matrix)
    tm_i=obj.shadow_transition_matrix(ii).Q;
    handle_struct.transition_matrix(ii,1)=rise_anonymous(tm_i.size,str2func(['@(y,x,ss,param)[',tm_i.string,']']),tm_i.indices);
end
%% steady state
% the steady state has equalities and will need to be evaluated. It can't
% be transformed into an anonymous function.
if isempty(obj.steady_state_shadow_model)
    handle_struct.steady_state_model=[];
else
    tmp=substitute_definitions(cell2mat(obj.steady_state_shadow_model'),defcell);
    handle_struct.steady_state_model=tmp;
    if obj.is_imposed_steady_state
        obj.is_stationary_model=true;
    end
end

%% model derivatives
JESP= obj.model_derivatives.Endogenous_Shocks_parameters;
aux_jac=cell2mat(obj.model_derivatives.auxiliary_jacobian(:)');
if isempty(aux_jac)
    aux_jac='';
end
JESP=rise_anonymous(JESP.size,str2func(['@(y,x,ss,param)[',JESP.string,']']),JESP.indices,aux_jac);
handle_struct.derivatives=JESP;


%% planner objective
% in this case there are items that could be written as anonymous
% functions, such as the objective, commiment and discount. We could write
% a function such as @(y,x,ss,param)deal(objective,commitment,discount),
% which would automatically return 3 outputs and would error if a number of
% outputs different from 3 was called. This drawback would be circumvented
% by writing to the toolbox an additional function calling the anonymous
% function as input and explicitly returning the outputs. We do not do this
% as there are also the derivatives which need to be evaluated. The
% solution is just to evaluate a long sequence of items and have a policy
% evaluation function as in old days.
if isempty(obj.planner.shadow_model) % function [objective,commitment,discount,der1,der2]=planner(y,ss,param)
    handle_struct.planner=[];
else
    items={'first_order_derivatives','second_order_derivatives','commitment','discount'};
    for ii=1:numel(items)
        vi=items{ii};
        tmp=obj.planner.(vi);
        handle_struct.planner.(vi)=rise_anonymous(tmp.size,str2func(['@(y,x,ss,param)[',tmp.string,']']),tmp.indices);
    end
end
obj.func_handles=handle_struct;
end

%% functions
function string=substitute_definitions(string,defcell)
for ii=1:size(defcell,1)
    string=strrep(string,defcell{ii,1},['(',defcell{ii,2},')']);
end
end

