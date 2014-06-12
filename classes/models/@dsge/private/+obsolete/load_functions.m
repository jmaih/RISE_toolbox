function obj=load_functions(obj)
if isempty(obj)
    obj=struct(); 
    return
end

MainFolder=obj.options.results_folder;
SubFoldersList={'graphs','estimation','simulations'}; %,'routines'

if ~exist(MainFolder,'dir')
    mkdir(MainFolder)
end
obj.folders_paths=struct();
for ifold=1:numel(SubFoldersList)
	subfolder=[MainFolder,filesep,SubFoldersList{ifold}];
	obj.folders_paths.(SubFoldersList{ifold})=subfolder;
	if ~exist(subfolder,'dir')
	    mkdir(subfolder)
	end
end

% 2- I HAVE CHANGE THE LEAD-LAG-INCIDENCE AND THIS MOST LIKELY HAS CONSEQUENCES

% likelihood functions
%---------------------

if obj.is_dsge_var_model
    obj.routines.likelihood=@likelihood_dsge_var;
elseif obj.is_optimal_simple_rule_model
    obj.routines.likelihood=@likelihood_optimal_simple_rule;
else
    obj.routines.likelihood=@likelihood_markov_switching_dsge;
end

%% steady state
% the steady state has equalities and will need to be evaluated. It can't
% be transformed into an anonymous function.
keyboard
if isempty(obj.steady_state_shadow_model)
    handle_struct.steady_state_model=[];
else
    steady_state_model=...
        substitute_definitions(cell2mat(obj.steady_state_shadow_model'),defcell);
    routines.steady_state_model=recreate_evaluation_structure(steady_state_model,obj.input_list,{'y','param'});
    if obj.is_imposed_steady_state
        obj.is_stationary_model=true;
    end
end

obj.func_handles=handle_struct;

    function batch=recreate_evaluation_structure(batch,input_list,ouput_list)
        if iscell(batch)
            batch=cell2mat(batch(:)');
        end
        if ~isstruct(batch)
            batch=struct('code',batch,'argins',{input_list},'argouts',{ouput_list});
        end
    end
end

%% functions

% function string=vectorize_string(string,varlist)
% % vectorized version
% if ischar(varlist)
%     varlist=cellstr(varlist);
% end
% varlist=varlist(:)';
% varlist=cell2mat(strcat(varlist,'|'));
% varlist=varlist(1:end-1);
% string=regexprep(string,'(?<!\.)([\^/*]{1})','.$1');
% string=regexprep(string,['(?<!\w)(',varlist,')(\()(\d+)(\))'],'$1$2$3,:$4');
% end

% function fhandle=create_function_handle(input_list,cellarray)
% if ischar(input_list)
%     input_list=cellstr(input_list);
% end
% input_list=input_list(:)';
% input_list=strcat(input_list,',');
% input_list=cell2mat(input_list);
% input_list=input_list(1:end-1);
% tmp=cell2mat(cellarray);
% tmp(isspace(tmp))=[];
% tmp(end)=[]; % remove last semicolon
% tmp=['@(',input_list,')[',tmp,']'];
% fhandle=str2func(tmp);
% end

% % dynamic: endogenous and exogenous
% %----------------------------------
% shd=obj.equations.shadow_dynamic;
% % replace all x's
% iter=nnz(obj.lead_lag_incidence);
% for ii=1:exo_nbr
%     iter=iter+1;
%     shd=strrep(shd,['x(',sprintf('%0.0f',ii),')'],['y(',sprintf('%0.0f',iter),')']);
% end
% % replace y by z
% shd=regexprep(shd,'(?<!\w)y(','z(');
% handle_struct.dynamic=utils.code.code2func(shd,{'z','x','ss','param','def','s0','s1'},true);
% 
% % vectorized form: endogenous and exogenous
% %------------------------------------------
% string=vectorize_string(shd,{'z'});
% tmp=utils.code.code2func(string,{'z','x','ss','param','def','s0','s1'},false);
% % since it has been vectorized, we need to modify the information about the
% % size and the number non-zero terms
% tmp.size(2)=nnz(obj.lead_lag_incidence)+exo_nbr;
% tmp.nnz_derivs=nan; % we don't know the number of non-zero terms
% handle_struct.vectorized_dynamic=tmp;
% 
% % parameters: non-vectorized form 
% %--------------------------------
% % the definitions have to be substituted unfortunately we are taking
% % derivatives wrt the parameters 
% shd=substitute_definitions(shd,defcell);
% handle_struct.dynamic_params=utils.code.code2func(shd,{'z','x','ss','param','def','s0','s1'},true);

% % parameters: vectorized form 
% %----------------------------
% % vectorization has to be done also with respect to z because some
% % equations may not have parameters and this would kill then concatenation.
% tmp=utils.code.code2func(string,{'z','x','ss','param','def','s0','s1'},false);
% % since it has been vectorized, we need to modify the information about the
% % size and the number non-zero terms
% tmp.size(2)=param_nbr;
% tmp.nnz_derivs=nan; % we don't know the number of non-zero terms
% handle_struct.vectorized_dynamic_params=tmp;

% %% planner objective
% if isempty(obj.planner_system.shadow_model) % function [objective,commitment,discount,der1,der2]=planner(y,ss,param)
%     handle_struct.planner=[];
% else
%     % this element is to be evaluated using the online function evaluator.
%     % it returns 5 outputs in the following order: loss, commitment,
%     % discount, hessian and jacobian
%     handle_struct.planner=obj.planner_system.LossComDiscHessJac;
% end
% 
% % handle_struct=set_structure_to_hard_function(handle_struct);
% % obj.model_derivatives=set_structure_to_hard_function(obj.model_derivatives);

