function obj=load_functions(obj)
if isempty(obj)
    obj=struct('rise_functions2disk',false);
    return
end

MainFolder=obj.options.results_folder;
SubFoldersList={'graphs','estimation','simulations','routines'};

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

% Get rid of definitions
defcell=substitute_definitions_in_definitions(obj.definitions);

% initialize
handle_struct=struct();

% function handles
if obj.is_dsge_var_model
    handle_struct.likelihood=@likelihood_dsge_var;
elseif obj.is_optimal_simple_rule_model
    handle_struct.likelihood=@likelihood_optimal_simple_rule;
else
    handle_struct.likelihood=@likelihood_markov_switching_dsge;
end
%% number of elements
exo_nbr=sum(obj.exogenous.number);
param_nbr=sum(obj.parameters.number);
eqtn_nbr=obj.equations.number;
%% anonymous function for definitions
theDef=obj.definitions.shadow_dynamic;
theDef=[{['def=zeros(',sprintf('%0.0f',obj.definitions.number),',1);']},theDef(:)'];
recreate_function_handle(theDef,{'param'},{'def'},'definitions');

%% static
mystatic=strcat('RES(',int2str((1:eqtn_nbr)'),')=',...
    obj.equations.shadow_static);
mystatic=[['RES=zeros(',sprintf('%0.0f',eqtn_nbr),',1);'];mystatic];
recreate_function_handle(mystatic,obj.input_list,{'RES'},'static');
% handle_struct.static=create_function_handle(obj.input_list,{obj.equations.shadow_static});

%% balanced growth path
mybalanced_growth=strcat('RES(',int2str((1:2*eqtn_nbr)'),')=',...
    obj.equations.shadow_balanced_growth_path);
mybalanced_growth=[['RES=zeros(',sprintf('%0.0f',2*eqtn_nbr),',1);'];mybalanced_growth];
recreate_function_handle(mybalanced_growth,obj.input_list,{'RES'},'balanced_growth');
% handle_struct.balanced_growth=create_function_handle(obj.input_list,tmp);

%% dynamic and vectorized dynamic
shd=obj.equations.shadow_dynamic;
% replace all x's
iter=nnz(obj.Lead_lag_incidence);
for ii=1:exo_nbr
    iter=iter+1;
    shd=strrep(shd,['x(',sprintf('%0.0f',ii),')'],['y(',sprintf('%0.0f',iter),')']);
end
% replace y by z
shd=regexprep(shd,'(?<!\w)y(','z(');

% endogenous and exogenous
% non-vectorized form
tmp=strcat('RES(',int2str((1:eqtn_nbr)'),')=',shd);
tmp=[['RES=zeros(',sprintf('%0.0f',eqtn_nbr),',1);'];tmp];
recreate_function_handle(tmp,{'z','x','ss','param','def','s0','s1'},{'RES'},'dynamic');

% vectorized form
string=vectorize_string(shd,{'z'});
tmp=strcat('RES(',int2str((1:eqtn_nbr)'),',:)=',string);
tmp=[['RES=zeros(',sprintf('%0.0f',eqtn_nbr),...
    ',',sprintf('%0.0f',nnz(obj.Lead_lag_incidence)+exo_nbr),');'];tmp];
recreate_function_handle(tmp,{'z','x','ss','param','def','s0','s1'},{'RES'},'vectorized_dynamic');

% for the parameters, the definitions have to be substituted
% unfortunately we are taking derivatives wrt the parameters
shd=substitute_definitions(shd,defcell);
% non-vectorized form
tmp=strcat('RES(',int2str((1:eqtn_nbr)'),')=',shd);
tmp=[['RES=zeros(',sprintf('%0.0f',eqtn_nbr),',1);'];tmp];
recreate_function_handle(tmp,{'z','x','ss','param','def','s0','s1'},{'RES'},'dynamic_params');

% vectorized form: NB: vectorization has to be done also with respect to z
% because some equations may not have parameters and this would kill the
% concatenation.
string=vectorize_string(shd,{'param','z'});
tmp=strcat('RES(',int2str((1:eqtn_nbr)'),',:)=',string);
tmp=[['RES=zeros(',sprintf('%0.0f',eqtn_nbr),...
    ',',sprintf('%0.0f',param_nbr),');'];tmp];
recreate_function_handle(tmp,{'z','x','ss','param','def','s0','s1'},{'RES'},'vectorized_dynamic_params');

%% transition matrix
kode=obj.shadow_transition_matrix;
argout=obj.shadow_transition_matrix.argouts{1};
kode.code=['retcode=0;',kode.code,'if any(isnan(',argout,'(:))) || any(',argout,'(:)<0) || any(',argout,'(:)>1);',argout,'=[];retcode=3;end;'];
% augment the code with expections. This could also have been done at
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

handle_struct=set_structure_to_hard_function(handle_struct);
obj.model_derivatives=set_structure_to_hard_function(obj.model_derivatives);

obj.func_handles=handle_struct;

    function handle_struct=set_structure_to_hard_function(handle_struct)
        this_folder=pwd;
        cd(obj.folders_paths.routines)
        myfields=fieldnames(handle_struct);
        for ifield=1:numel(myfields)
            cell_flag=iscell(handle_struct.(myfields{ifield}));
            if ~((cell_flag && isstruct(handle_struct.(myfields{ifield}){1}))||...
                    (~cell_flag && isstruct(handle_struct.(myfields{ifield}))))
                continue
            end
            this_name=[myfields{ifield},'___'];
            if cell_flag
                code2file(handle_struct.(myfields{ifield}){1},this_name);
                handle_struct.(myfields{ifield}){1}=str2func(['@',this_name]);
            else
                code2file(handle_struct.(myfields{ifield}),this_name);
                handle_struct.(myfields{ifield})=str2func(['@',this_name]);
            end
        end
        cd(this_folder)
        builtin('rehash');
    end

    function recreate_function_handle(batch,input_list,ouput_list,fname)
        if iscell(batch)
            batch=cell2mat(batch(:)');
        end
        if ~isstruct(batch)
            batch=struct('code',batch,'argins',{input_list},'argouts',{ouput_list});
        end
        if obj.options.rise_functions2disk
            routine_dir_name=[obj.options.results_folder,filesep,'routines'];
            code2file(batch,fname)
            movefile([fname,'.m'],routine_dir_name)
            curr_dir=pwd;
            % go into the folder and get the handle
            cd(routine_dir_name);
            handle_struct.(fname)=str2func(['@',fname]);
            cd(curr_dir);
        else
            handle_struct.(fname)=batch;
        end
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
