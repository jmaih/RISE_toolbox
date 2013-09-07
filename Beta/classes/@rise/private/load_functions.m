function obj=load_functions(obj)
defaultOptions=struct('rise_functions2disk',false);
if isempty(obj)
    obj=defaultOptions;
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
if obj.is_svar_model
    return
elseif obj.is_dsge_var_model
    handle_struct.likelihood=@likelihood_dsge_var;
elseif obj.is_optimal_simple_rule_model
    handle_struct.likelihood=@likelihood_optimal_simple_rule;
else
    handle_struct.likelihood=@likelihood_markov_switching_dsge;
end

%% anonymous function for definitions
theDef={obj.definitions.shadow_dynamic};
theDef=[{['def=zeros(',int2str(numel(obj.definitions)),',1);']},theDef];
recreate_function_handle(theDef,{'param'},{'def'},'definitions');

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
mystatic=strcat('RES(',int2str((1:obj.NumberOfEquations)'),')=',...
    {obj.equations.shadow_static}');
mystatic=[['RES=zeros(',int2str(obj.NumberOfEquations),',1);'];mystatic];
recreate_function_handle(mystatic,obj.input_list,{'RES'},'static');
% handle_struct.static=create_function_handle(obj.input_list,{obj.equations.shadow_static});

%% balanced growth path
mybalanced_growth=cellstr(obj.equations(1).shadow_balanced_growth_path)';
for ii=2:obj.NumberOfEquations
    mybalanced_growth=[mybalanced_growth,cellstr(obj.equations(ii).shadow_balanced_growth_path)']; %#ok<AGROW>
end
mybalanced_growth=strcat('RES(',int2str((1:2*obj.NumberOfEquations)'),')=',...
    mybalanced_growth');
mybalanced_growth=[['RES=zeros(',int2str(2*obj.NumberOfEquations),',1);'];mybalanced_growth];
recreate_function_handle(mybalanced_growth,obj.input_list,{'RES'},'balanced_growth');
% handle_struct.balanced_growth=create_function_handle(obj.input_list,tmp);

%% dynamic and vectorized dynamic
shd={obj.equations.shadow_dynamic};
% replace all x's
iter=nnz(obj.Lead_lag_incidence);
for ii=1:obj.NumberOfExogenous
    iter=iter+1;
    shd=strrep(shd,['x(',int2str(ii),')'],['y(',int2str(iter),')']);
end
% replace y by z
shd=regexprep(shd,'(?<!\w)y(','z(');

% endogenous and exogenous
% non-vectorized form
tmp=strcat('RES(',int2str((1:obj.NumberOfEquations)'),')=',...
    shd');
tmp=[['RES=zeros(',int2str(obj.NumberOfEquations),',1);'];tmp];
recreate_function_handle(tmp,{'z','param','ss','def'},{'RES'},'dynamic');

% vectorized form
string=vectorize_string(shd,{'z'});
tmp=strcat('RES(',int2str((1:obj.NumberOfEquations)'),',:)=',...
    string');
tmp=[['RES=zeros(',int2str(obj.NumberOfEquations),...
    ',',int2str(nnz(obj.Lead_lag_incidence)+obj.NumberOfExogenous),');'];tmp];
recreate_function_handle(tmp,{'z','param','ss','def'},{'RES'},'vectorized_dynamic');

% for the parameters, the definitions have to be substituted
% unfortunately we are taking derivatives wrt the parameters
shd=substitute_definitions(shd,defcell);
% non-vectorized form
tmp=strcat('RES(',int2str((1:obj.NumberOfEquations)'),')=',...
    shd');
tmp=[['RES=zeros(',int2str(obj.NumberOfEquations),',1);'];tmp];
recreate_function_handle(tmp,{'param','z','ss'},{'RES'},'dynamic_params');

% vectorized form: NB: vectorization has to be done also with respect to z
% because some equations may not have parameters and this would kill the
% concatenation.
string=vectorize_string(shd,{'param','z'});
tmp=strcat('RES(',int2str((1:obj.NumberOfEquations)'),',:)=',...
    string');
tmp=[['RES=zeros(',int2str(obj.NumberOfEquations),...
    ',',int2str(obj.NumberOfParameters),');'];tmp];
recreate_function_handle(tmp,{'param','z','ss'},{'RES'},'vectorized_dynamic_params');

%% transition matrix
kode=obj.shadow_transition_matrix;
kode.code=['retcode=0;',kode.code,'if any(any(isnan(Q))) || any(any(Q<0)) || any(any(Q>1));Q=[];retcode=3;end;'];
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
    recreate_function_handle(steady_state_model,obj.input_list,{'y'},'steady_state_model');
    if obj.is_imposed_steady_state
        obj.is_stationary_model=true;
    end
end

%% model derivatives: 
% those elements are already in the adequate structure form
recreate_function_handle(obj.model_derivatives.Endogenous_Shocks,[],[],...
    'endo_exo_derivatives');

recreate_function_handle(obj.model_derivatives.StaticEndogenous,[],[],...
    'static_model_derivatives');

recreate_function_handle(obj.model_derivatives.Static_BGP_Endogenous,[],[],...
    'static_bgp_model_derivatives');

recreate_function_handle(obj.model_derivatives.Parameters,[],[],'param_derivatives');

%% planner objective
if isempty(obj.planner.shadow_model) % function [objective,commitment,discount,der1,der2]=planner(y,ss,param)
    handle_struct.planner=[];
else
    % this element is to be evaluated using the online function evaluator.
    % it returns 5 outputs in the following order: loss, commitment,
    % discount, hessian and jacobian
    handle_struct.planner=obj.planner.LossComDiscHessJac;
end
obj.func_handles=handle_struct;

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

function string=substitute_definitions(string,defcell)
for ii=1:size(defcell,1)
    string=strrep(string,defcell{ii,1},['(',defcell{ii,2},')']);
end
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
