function obj=load_functions(obj)
if isempty(obj)
    obj=struct(); %<-- obj=struct('rise_functions2disk',false);% this option is obsolete by now
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

% number of elements
%-------------------
exo_nbr=sum(obj.exogenous.number);
param_nbr=sum(obj.parameters.number);
% eqtn_nbr=obj.equations.number;

% Get rid of definitions: as they can be functions of themselves and hence
% do not lend themselves to be written as a function
defcell=substitute_definitions_in_definitions(obj.definitions);

% initialize handles
%-------------------
handle_struct=struct();

% add the definitions
%--------------------
handle_struct.definitions=code2func(defcell(:,2),'param',true);

% likelihood functions
%---------------------
if obj.is_dsge_var_model
    handle_struct.likelihood=@likelihood_dsge_var;
elseif obj.is_optimal_simple_rule_model
    handle_struct.likelihood=@likelihood_optimal_simple_rule;
else
    handle_struct.likelihood=@likelihood_markov_switching_dsge;
end

% static model
%-------------
handle_struct.static=code2func(obj.equations.shadow_static,obj.input_list,true);

% balanced growth path
%---------------------
handle_struct.balanced_growth=code2func(obj.equations.shadow_balanced_growth_path,obj.input_list,true);

% dynamic: endogenous and exogenous
%----------------------------------
shd=obj.equations.shadow_dynamic;
% replace all x's
iter=nnz(obj.Lead_lag_incidence);
for ii=1:exo_nbr
    iter=iter+1;
    shd=strrep(shd,['x(',sprintf('%0.0f',ii),')'],['y(',sprintf('%0.0f',iter),')']);
end
% replace y by z
shd=regexprep(shd,'(?<!\w)y(','z(');
handle_struct.dynamic=code2func(shd,{'z','x','ss','param','def','s0','s1'},true);

% vectorized form: endogenous and exogenous
%------------------------------------------
string=vectorize_string(shd,{'z'});
tmp=code2func(string,{'z','x','ss','param','def','s0','s1'},false);
% since it has been vectorized, we need to modify the information about the
% size and the number non-zero terms
tmp.size(2)=nnz(obj.Lead_lag_incidence)+exo_nbr;
tmp.nnz_derivs=nan; % we don't know the number of non-zero terms
handle_struct.vectorized_dynamic=tmp;

% parameters: non-vectorized form 
%--------------------------------
% the definitions have to be substituted unfortunately we are taking
% derivatives wrt the parameters 
shd=substitute_definitions(shd,defcell);
handle_struct.dynamic_params=code2func(shd,{'z','x','ss','param','def','s0','s1'},true);

% parameters: vectorized form 
%----------------------------
% vectorization has to be done also with respect to z because some
% equations may not have parameters and this would kill then concatenation.
tmp=code2func(string,{'z','x','ss','param','def','s0','s1'},false);
% since it has been vectorized, we need to modify the information about the
% size and the number non-zero terms
tmp.size(2)=param_nbr;
tmp.nnz_derivs=nan; % we don't know the number of non-zero terms
handle_struct.vectorized_dynamic_params=tmp;

%% transition matrix
kode=obj.shadow_transition_matrix;
argout=obj.shadow_transition_matrix.argouts{1};
kode.code=['retcode=0;',kode.code,'if any(isnan(',argout,'(:))) || any(',argout,'(:)<0) || any(',argout,'(:)>1);',argout,'=[];retcode=3;end;'];
% augment the code with exceptions. This could also have been done at
% compilation time
kode.argouts=[kode.argouts,'retcode'];
recreate_function_handle(kode,[],[],'transition_matrix');
% finally load it                           
%% steady state
% the steady state has equalities and will need to be evaluated. It can't
% be transformed into an anonymous function.
if isempty(obj.steady_state_shadow_model)
    handle_struct.steady_state_model=[];
else
    steady_state_model=...
        substitute_definitions(cell2mat(obj.steady_state_shadow_model'),defcell);
    recreate_function_handle(steady_state_model,obj.input_list,{'y','param'},'steady_state_model');
    if obj.is_imposed_steady_state
        obj.is_stationary_model=true;
    end
end

%% planner objective
if isempty(obj.planner_system.shadow_model) % function [objective,commitment,discount,der1,der2]=planner(y,ss,param)
    handle_struct.planner=[];
else
    % this element is to be evaluated using the online function evaluator.
    % it returns 5 outputs in the following order: loss, commitment,
    % discount, hessian and jacobian
    handle_struct.planner=obj.planner_system.LossComDiscHessJac;
end

% handle_struct=set_structure_to_hard_function(handle_struct);
% obj.model_derivatives=set_structure_to_hard_function(obj.model_derivatives);

obj.func_handles=handle_struct;

    function recreate_function_handle(batch,input_list,ouput_list,fname)
        if iscell(batch)
            batch=cell2mat(batch(:)');
        end
        if ~isstruct(batch)
            batch=struct('code',batch,'argins',{input_list},'argouts',{ouput_list});
        end
%         if obj.options.rise_functions2disk
%             routine_dir_name=[obj.options.results_folder,filesep,'routines'];
%             code2file(batch,fname)
%             movefile([fname,'.m'],routine_dir_name)
%             curr_dir=pwd;
%             % go into the folder and get the handle
%             cd(routine_dir_name);
%             handle_struct.(fname)=str2func(['@',fname]);
%             cd(curr_dir);
%         else
            handle_struct.(fname)=batch;
%         end
    end
end

%% functions

function string=vectorize_string(string,varlist)
% vectorized version
if ischar(varlist)
    varlist=cellstr(varlist);
end
varlist=varlist(:)';
varlist=cell2mat(strcat(varlist,'|'));
varlist=varlist(1:end-1);
string=regexprep(string,'(?<!\.)([\^/*]{1})','.$1');
string=regexprep(string,['(?<!\w)(',varlist,')(\()(\d+)(\))'],'$1$2$3,:$4');
end

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
