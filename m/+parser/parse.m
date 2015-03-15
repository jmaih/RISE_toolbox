function dictionary=parse(FileName,varargin)
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

% - can declare exogenous and make them observable at the
% same time. The exogenous that are observed are determisitic. This opens
% the door for estimating partial equilibrium models
% - the rise_flags are either a structure or a cell array with two
% columns!!! 

DefaultOptions=...
    struct('definitions_in_param_differentiation',false,...
    'rise_save_macro',false,...
    'max_deriv_order',2,'definitions_inserted',false,...
    'parse_debug',false);
DefaultOptions=utils.miscellaneous.mergestructures(DefaultOptions,...
    parser.preparse());
if nargin<1
    dictionary=DefaultOptions;
    return
end

%%
lv=length(varargin);
if rem(lv,2)~=0
    error('arguments must come in pairs')
end
if lv
    ccc=reshape(varargin(:),2,[]); clear varargin
    for ii=1:0.5*lv
        str=ccc{1,ii};
        if ~isfield(DefaultOptions,str)
            error([str,' Not an option for ',mfilename])
        end
        DefaultOptions.(str)=ccc{2,ii};
    end
end
%% general initializations

dictionary = parser.initialize_dictionary();
dictionary.definitions_inserted=DefaultOptions.definitions_inserted;
dictionary.parse_debug=DefaultOptions.parse_debug;
%% set various blocks

% first output: the dictionary.filename

FileName(isspace(FileName))=[];
loc=strfind(FileName,'.');
if isempty(loc)
    dictionary.filename=FileName;
    if exist([FileName,'.rs'],'file')
        FileName=[FileName,'.rs'];
    elseif exist([FileName,'.rz'],'file')
        FileName=[FileName,'.rz'];
    elseif exist([FileName,'.dsge'],'file')
        FileName=[FileName,'.dsge'];
    else
        error([mfilename,':: ',FileName,'.rs or ',FileName,'.rz  or ',FileName,'.dsge not found'])
    end
else
    ext=FileName(loc:end);
    %     if ~(strcmp(ext,'.rs')||strcmp(ext,'.dyn')||strcmp(ext,'.mod'))
    if ~ismember(ext,{'.rs','.rz','.dsge'})%||strcmp(ext,'.dyn')||strcmp(ext,'.mod'))
        error([mfilename,':: Input file is expected to be a .rs or a .rz  or a .dsge file'])
    end
    dictionary.filename=FileName(1:loc-1);
end

% read file and remove comments
% RawFile=read_file(FileName,DefaultOptions.rise_flags);
[RawFile,has_macro]=parser.preparse(FileName,DefaultOptions.rise_flags);

% write the expanded version
if has_macro && DefaultOptions.rise_save_macro
    newfile='';
    thedot=strfind(FileName,'.');
    fid=fopen([FileName(1:thedot-1),'_expanded.dsge'],'w');
    for irow=1:size(RawFile,1)
        write_newfilename=isempty(newfile)||~strcmp(newfile,RawFile{irow,2});
        if write_newfilename
            newfile=RawFile{irow,2};
            fprintf(fid,'%s\n',['// ',newfile,' line ',sprintf('%0.0f',RawFile{irow,3})]);
        end
        fprintf(fid,'%s\n',RawFile{irow,1});
    end
    fclose(fid);
end

[blocks,dictionary.markov_chains]=parser.file2blocks(RawFile);

dictionary.raw_file=cell2struct(RawFile,{'code','filename','line_numbers'},2);
dictionary.rawfile_triggers={blocks.trigger};

clear RawFile

%% Populate the dictionary
[dictionary,blocks]=parser.declarations2dictionary(dictionary,blocks);


%% add the chain names to dictionary so that the parsing of the model does
% not complain although the chain names acquired so far through the
% parameters are not expected to enter the model block
dictionary.chain_names={dictionary.markov_chains.name};

%% turn all log_vars into exp(log_vars) both in the model block
% and update the names of the variables in the process
[dictionary,blocks,old_endo_names]=parser.logvars2logvars(dictionary,blocks);

%% Model block
% now with the endogenous, exogenous, parameters in hand, we can process
current_block_id=find(strcmp('model',{blocks.name}));
more_string='';
if dictionary.definitions_inserted
    more_string='(& possibly definitions insertions)';
end
if dictionary.parse_debug
    profile off
    profile on
else
    tic
end
[Model_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'model');
if dictionary.parse_debug
    profile off
    profile viewer
    keyboard
else
    disp([mfilename,':: Model parsing ',more_string,'. ',sprintf('%0.4f',toc),' seconds'])
end

if isempty(Model_block)
    error([mfilename,':: no model declared'])
end
% remove item from block
blocks(current_block_id)=[];
%% after parsing the model block, update the markov chains (time-varying probabilities)
% sorting the endogenous switching probabilities is more or less useless
dictionary.time_varying_probabilities=sort(dictionary.time_varying_probabilities);

%% add more equations if necessary, such that the max_lag and max_lead are 1
% replace exogenous lags with auxiliary endogenous variables
% update the number of endogenous variables accordingly
% Also keep the same equations to be added to the steady state model
auxiliary_steady_state_equations=cell(0,4);
orig_endogenous_current=dictionary.orig_endogenous; % these are the variables without the augmentation add-on
all_fields=fieldnames(parser.listing);
main_endo_fields={'name','tex_name','max_lead','max_lag','is_log_var',...
    'is_auxiliary','is_trans_prob'};
useless_fields=setdiff(all_fields,main_endo_fields);
overall_max_lead_lag=max([dictionary.orig_endogenous.max_lead]);
overall_max_lead_lag=max([abs([dictionary.orig_endogenous.max_lag]),overall_max_lead_lag]);
for ii=1:numel(dictionary.orig_endogenous)
    if dictionary.orig_endogenous(ii).max_lead<2
        continue
    end
    vname=dictionary.orig_endogenous(ii).name;
    lead_i=dictionary.orig_endogenous(ii).max_lead;
    vold=vname;
    for i2=2:lead_i
        new_name=parser.create_auxiliary_name(vname,i2-1);
        new_var=struct('name',new_name,'tex_name','','max_lead',1,...
            'max_lag',0,'is_log_var',false,...
            'is_auxiliary',true,'is_trans_prob',false);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vold,1}',{';',[]}'],0,1,'normal'}]; %#ok<*AGROW> %
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{vold,0}',{';',[]}'],0,0,'normal'}]; %<-- add maxLag and maxLead
        if i2-1==1 % set the lead to 1
            dictionary.orig_endogenous(ii).max_lead=1;
        end
        for ifield=1:numel(useless_fields)
            new_var.(useless_fields{ifield})=nan;
        end
        dictionary.orig_endogenous=[dictionary.orig_endogenous,new_var];
        orig_endogenous_current=[orig_endogenous_current,dictionary.orig_endogenous(ii)];
        vold=new_var.name;
    end
end

% now merge endogenous and exogenous
n_endog=numel(dictionary.orig_endogenous);
variables=[dictionary.orig_endogenous,dictionary.exogenous];
for ii=1:numel(variables)
    is_endogenous=ii<=n_endog;
    % for endogenous,lags have to be greater than 1, for exogenous, lags
    % have to be greater than 0
    if (is_endogenous && variables(ii).max_lag>-2)||(~is_endogenous && variables(ii).max_lag>-1)
        continue
    end
    vname=variables(ii).name;
    if ~is_endogenous
        % then create an auxiliary variable and an auxiliary equation
        % before proceeding
        new_var=struct('name',[vname,'_0'],'tex_name','','max_lead',0,...
            'max_lag',-1,'is_log_var',false,'is_auxiliary',true,...
            'is_trans_prob',false);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vname,0}',{';',[]}'],0,0,'normal'}]; %
        % The steady state is computed with zero shocks. and so instead of
        % setting the name of the shock, I write 0
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{'0',0}',{';',[]}'],0,0,'normal'}]; %<-- add maxLag and maxLead
        for ifield=1:numel(useless_fields)
            new_var.(useless_fields{ifield})=nan;
        end
        dictionary.orig_endogenous=[dictionary.orig_endogenous,new_var];
        orig_endogenous_current=[orig_endogenous_current,variables(ii)];
        % correct the lag structure of the original variable
        vloc=strcmp(variables(ii).name,{dictionary.exogenous.name});
        dictionary.exogenous(vloc).max_lag=0;
        vname=new_var.name;
    end
    vold=vname;
    for i2=2:abs(variables(ii).max_lag)
        % update the lag structure of the original variable if it is
        % endogenous
        if i2==2 && is_endogenous
            vloc=strcmp(variables(ii).name,{dictionary.orig_endogenous.name});
            dictionary.orig_endogenous(vloc).max_lag=-1;
        end
        new_name=parser.create_auxiliary_name(vname,uminus(i2-1));
        new_var=struct('name',new_name,'tex_name','','max_lead',0,...
            'max_lag',-1,'is_log_var',false,'is_auxiliary',true,...
            'is_trans_prob',false);
        
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vold,-1}',{';',[]}'],-1,0,'normal'}]; %
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{vold,0}',{';',[]}'],0,0,'normal'}]; %<-- add maxLag and maxLead
        for ifield=1:numel(useless_fields)
            new_var.(useless_fields{ifield})=nan;
        end
        dictionary.orig_endogenous=[dictionary.orig_endogenous,new_var];
        orig_endogenous_current=[orig_endogenous_current,variables(ii)];
        vold=new_var.name;
    end
end
% remove this from workspace
clear variables

%% Now we can re-order the original (and augmented) endogenous
[~,tags]=sort({dictionary.orig_endogenous.name});
dictionary.orig_endogenous=dictionary.orig_endogenous(tags);
% as well as the same list but without the augmentation add-ons
orig_endogenous_current=orig_endogenous_current(tags);

%% Now we re-write the model and update leads and lags
% at the same time, construct the incidence and occurrence matrices
orig_endo_nbr=numel(dictionary.orig_endogenous);

number_of_equations=size(Model_block,1);
Occurrence=false(number_of_equations,orig_endo_nbr,3);
equation_type=ones(number_of_equations,1);
for ii=1:number_of_equations
    eq_i= Model_block{ii,1};
    maxLag= Model_block{ii,2};
    maxLead= Model_block{ii,3};
    if strcmp(Model_block{ii,4},'def') % <---ismember(eq_i{1,1},dictionary.definitions) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=2;
    elseif strcmp(Model_block{ii,4},'tvp') % <---ismember(eq_i{1,1},dictionary.time_varying_probabilities) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=3;
    elseif strcmp(Model_block{ii,4},'mcp') % <---ismember(eq_i{1,1},dictionary.time_varying_probabilities) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=4;
    end
    for i2=1:size(eq_i,2)
        if ~isempty(eq_i{2,i2})
            vname=eq_i{1,i2};
            status=dictionary.determine_status(vname,dictionary);
            time=-1; new_item=false;
            if ~isempty(eq_i{2,i2})&& abs(eq_i{2,i2})>0
                if strcmp(status,'param')
                    if eq_i{2,i2}>1||eq_i{2,i2}<0
                        error('parameters can only have leads of 1')
                    end
                elseif strcmp(status,'x')
                    new_item=true;
                    if abs(eq_i{2,i2})==1
                        % change the name, but not the lag structure
                        eq_i{1,i2}=[vname,'_0'];
                    else
                        % change both the name and the lag structure
                        eq_i{1,i2}=parser.create_auxiliary_name([vname,'_0'],eq_i{2,i2}+1);
                    end
                elseif strcmp(status,'y')&& abs(eq_i{2,i2})>1
                    new_item=true;
                    % change both the name and the lag structure
                    if eq_i{2,i2}>0
                        eq_i{1,i2}=parser.create_auxiliary_name(vname,eq_i{2,i2}-1);
                        time=1;
                    elseif eq_i{2,i2}<0
                        eq_i{1,i2}=parser.create_auxiliary_name(vname,eq_i{2,i2}+1);
                    end
                end
            end
            if new_item
                eq_i{2,i2}=time;
            end
            var_loc=strcmp(eq_i{1,i2},{dictionary.orig_endogenous.name});
            lag_or_lead=eq_i{2,i2}+2;
            if any(var_loc) && equation_type(ii)==2
                error([mfilename,':: equation (',sprintf('%0.0f',ii),') detected to be a definition cannot contain variables'])
            elseif equation_type(ii)==3 && ismember(lag_or_lead,[3,1])
                error([mfilename,':: equation (',sprintf('%0.0f',ii),') detected to describe endogenous switching cannot contain leads or lags'])
            end
            Occurrence(ii,var_loc,lag_or_lead)=true;
        end
    end
    Model_block{ii,1}= eq_i;
    Model_block{ii,2}= max(-1,maxLag);
    Model_block{ii,3}= min(1,maxLead);
end
% keep only the structural equations
Occurrence=Occurrence(equation_type==1,:,:);

%% Steady state Model block
% now with the endogenous, exogenous, parameters in hand, we can process
% the steady state model block

static=struct('is_imposed_steady_state',false,'is_unique_steady_state',false);
current_block_id=find(strcmp('steady_state_model',{blocks.name}));
if dictionary.parse_debug
    profile off
    profile on
else
    tic
end
[SteadyStateModel_block,dictionary,static]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'steady_state_model',static);
if dictionary.parse_debug
    profile off
    profile viewer
    keyboard
else
    disp([mfilename,':: Steady State Model parsing . ',sprintf('%0.4f',toc),' seconds'])
end

% remove item from block
blocks(current_block_id)=[];

% if there are steady state equations, then add the auxiliary equations
if ~isempty(SteadyStateModel_block)
    SteadyStateModel_block=[SteadyStateModel_block
        auxiliary_steady_state_equations];
end

%% exogenous definitions block
% this block needs a special treatment as the information provided here
% will only be used when building the dataset for estimation and/or during
% forecasting.
current_block_id=find(strcmp('exogenous_definition',{blocks.name}));
if dictionary.parse_debug
    profile off
    profile on
else
    tic
end
[ExogenousDefinition_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'exogenous_definition');
if dictionary.parse_debug
    profile off
    profile viewer
    keyboard
else
    disp([mfilename,':: exogenous definitions block parsing . ',sprintf('%0.4f',toc),' seconds'])
end

% remove item from block
blocks(current_block_id)=[];
% the equations have been validated, now rebuild them and keep a list of
% the variables defined
DefinedExoList=cell(1,0);%{}
for ii=1:size(ExogenousDefinition_block,1)
    DefinedExoList=[DefinedExoList,ExogenousDefinition_block{ii,1}(1,1)];
    eq_i='';
    for jj=1:size(ExogenousDefinition_block{ii,1},2)
        eq_i=[eq_i,ExogenousDefinition_block{ii,1}{1,jj}];
        if ~isempty(ExogenousDefinition_block{ii,1}{2,jj}) && ExogenousDefinition_block{ii}{2,jj}~=0
            eq_i=[eq_i,'{',sprintf('%0.0f',ExogenousDefinition_block{ii,1}{2,jj}),'}'];
        end
    end
    ExogenousDefinition_block{ii,1}=eq_i;
end
% assign information to dictionary
dictionary.exogenous_equations=struct('name',DefinedExoList,'equation',transpose(ExogenousDefinition_block(:,1)));
clear ExogenousDefinition_block DefinedExoList
%% optimal policy block
current_block_id=find(strcmp('planner_objective',{blocks.name}));
[dictionary,PlannerObjective_block,is_model_with_planner_objective]=...
    parser.planner_objective(dictionary,blocks(current_block_id).listing);

% remove item from block
blocks(current_block_id)=[];

%% parameterization block
current_block_id=find(strcmp('parameterization',{blocks.name}));

dictionary.Parameterization_block=parser.capture_parameterization(dictionary,blocks(current_block_id).listing);

% remove item from block
blocks(current_block_id)=[];

%% parameter restrictions block.

current_block_id=find(strcmp('parameter_restrictions',{blocks.name}));

[Param_rest_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'parameter_restrictions');
% remove item from block
blocks(current_block_id)=[]; %#ok<NASGU>

% remove the columns with information about the maxLead and maxLag: the
% capture of parameterization does not require it and might even crash
dictionary.Param_rest_block=Param_rest_block(:,1);
clear Param_rest_block

%% Lump together the model and steady-state model
AllModels=[Model_block;
    SteadyStateModel_block;
    PlannerObjective_block];
ss_eq_nbr=size(SteadyStateModel_block,1);
% steady state equations are identified by number 5
planobj_eq_nbr=size(PlannerObjective_block,1);
% dictionary.planner_system objective equations are identified by number 6
equation_type=[equation_type;5*ones(ss_eq_nbr,1);6*ones(planobj_eq_nbr,1)];

% collect information about leads and lags which will be used to determine
% the status of the lagrange multipliers later on. But keep only the
% information about the structural equations
%--------------------------------------------------------------------------
equations_maxLag_maxLead=cell2mat(Model_block(equation_type==1,2:3));

clear Model_block SteadyStateModel_block PlannerObjective_block
%% Incidence matrix
% Now and only now can we build the incidence matrix
dictionary.lead_lag_incidence.before_solve=zeros(3,orig_endo_nbr);
for ii=1:3
    for jj=1:orig_endo_nbr
        if any(Occurrence(:,jj,ii))
            dictionary.lead_lag_incidence.before_solve(ii,jj)=1;
        end
    end
end
dictionary.lead_lag_incidence.before_solve=transpose(flipud(dictionary.lead_lag_incidence.before_solve));
dictionary.lead_lag_incidence.before_solve(dictionary.lead_lag_incidence.before_solve>0)=1:nnz(dictionary.lead_lag_incidence.before_solve);

appear_as_current=dictionary.lead_lag_incidence.before_solve(:,2)>0;
if any(~appear_as_current)
    disp('The following variables::')
    disp(old_endo_names(~appear_as_current))
    error('do not appear as current')
end

%% models in shadow/technical/tactical form

% extract the mcp equations as they have to be processed after the final
% list of endogenous is known...

is_mcp=equation_type==4;
mcp_eqtns=AllModels(is_mcp,:);
mcp_type=equation_type(is_mcp);
equation_type(is_mcp)=[];
AllModels(is_mcp,:)=[];

[dictionary,...
    dynamic,...
    stat,...
    defs,...
    shadow_tvp]=parser.shadowize(dictionary,AllModels,...
    equation_type,orig_endogenous_current,overall_max_lead_lag);
orig_definitions=defs.original;
shadow_definitions=defs.shadow;
static=utils.miscellaneous.mergestructures(static,stat);

% collect the future switching parameters
%----------------------------------------
fsp=regexp(dynamic.shadow_model,'sparam\(\d+\)','match');
fsp=unique([fsp{:}]);
switching_parameters_leads_index=regexprep(fsp,'sparam\((\d+)\)','$1');
switching_parameters_leads_index=str2num(char(switching_parameters_leads_index)); %#ok<ST2NM>
dictionary.fsp=fsp;

static.steady_state_shadow_model=strrep(static.steady_state_shadow_model,...
    'x1_=','[x1_,fval,exitflag]=');
static.steady_state_shadow_model=strrep(static.steady_state_shadow_model,...
    ',x0_,',',x0_,options,');
old_shadow_steady_state_model=static.steady_state_shadow_model;
static.steady_state_shadow_model=cell(0,1);
fsolve_nbr=0;
for ii=1:numel(old_shadow_steady_state_model)
    eq_i=old_shadow_steady_state_model{ii};
    if ~isempty(strfind(eq_i,'argzero'))
        eq_i=strrep(eq_i,'argzero','fsolve');
        static.steady_state_shadow_model=[static.steady_state_shadow_model;{eq_i}];
        eq_i={'retcode=1-(exitflag==1);'};
        if ii<numel(old_shadow_steady_state_model)
            eq_i=[eq_i;{'if ~retcode,'}];
            fsolve_nbr=fsolve_nbr+1;
        end
        static.steady_state_shadow_model=[static.steady_state_shadow_model;eq_i];
    else
        static.steady_state_shadow_model=[static.steady_state_shadow_model;{eq_i}];
    end
end
for ii=1:fsolve_nbr
    static.steady_state_shadow_model=[static.steady_state_shadow_model;{'end;'}];
end
clear old_shadow_steady_state_model
% now add the prelude
if ~isempty(static.steady_state_shadow_model)
    static.steady_state_shadow_model=[{['y=zeros(',sprintf('%0.0f',orig_endo_nbr),',1);']};static.steady_state_shadow_model];
end
%% replace the list of definition names with definition equations
dictionary.definitions=struct('model',orig_definitions,'shadow',shadow_definitions);
%% initialize the routines structure
routines=struct();

%% definitions routine
routines.definitions=utils.code.code2func(...
    parser.substitute_definitions(shadow_definitions),'param');
%% steady state model routine (cannot be written as a function)
routines.steady_state_model=parser.substitute_definitions(...
    static.steady_state_shadow_model,shadow_definitions);
routines.steady_state_model=struct('code',cell2mat(routines.steady_state_model(:)'),...
    'argins',{parser.input_list},...
    'argouts',{{'y','param'}});
%% load transition matrices and transform markov chains
% the Trans mat will go into the computation of derivatives
% dictionary=parser.transition_probabilities(dictionary,shadow_tvp);
bug_fix_shadow_defs=[];
if ~isempty(dictionary.definitions)
    bug_fix_shadow_defs=dictionary.definitions.shadow;
end
probability_of_commitment=[];
if ~isempty(dictionary.planner_system.shadow_model)
    probability_of_commitment=dictionary.planner_system.shadow_model{2};
    probability_of_commitment=strrep(probability_of_commitment,'commitment-','');
    probability_of_commitment=strrep(probability_of_commitment,';','');
    probability_of_commitment=probability_of_commitment(2:end-1);
end
[~,...
    routines.transition_matrix,...
    dictionary.markov_chains,...
    myifelseif...
    ]=parser.transition_probabilities(...
    dictionary.input_list,...
    {dictionary.parameters.name},dictionary.markov_chains,...
    shadow_tvp,bug_fix_shadow_defs,probability_of_commitment);
% flag for endogenous probabilities
%---------------------------------
dictionary.is_endogenous_switching_model=any(dictionary.markov_chains.chain_is_endogenous);
%% Computation of derivatives
%
disp(' ')
disp('Now computing symbolic derivatives...')
disp(' ')
max_deriv_order=max(1,DefaultOptions.max_deriv_order);
exo_nbr=numel(dictionary.exogenous);
% switching_ones=find([dictionary.parameters.is_switching]);

% dynamic model wrt endogenous, exogenous and leads of switching parameters
%--------------------------------------------------------------------------
% no leads or lags in the endogenous probabilities: then we know that the
% status (static,pred,both,frwrd) of the variables does not change.

% v=[f+,b+,s0,p_0,b_0,f_0,p_minus,b_minus,e_0,theta_plus]
routines.dynamic=utils.code.code2func(dynamic.shadow_model);

[wrt,dictionary.v,...
    dictionary.locations.before_solve,...
    dictionary.siz.before_solve,...
    dictionary.order_var.before_solve,...
    dictionary.inv_order_var.before_solve,...
    dictionary.steady_state_index]=dynamic_differentiation_list(...
    dictionary.lead_lag_incidence.before_solve,exo_nbr,switching_parameters_leads_index);
%----------------------
if dictionary.parse_debug
    profile off
    profile on
end
routines.probs_times_dynamic=parser.burry_probabilities(dynamic.shadow_model,myifelseif);
[routines.probs_times_dynamic_derivatives,numEqtns,numVars,jac_toc,...
    original_funcs]=differentiate_system(...
    routines.probs_times_dynamic,...
    dictionary.input_list,...
    wrt,...
    max_deriv_order);
routines.probs_times_dynamic=utils.code.code2func(routines.probs_times_dynamic);
if dictionary.parse_debug
    profile off
    profile viewer
    keyboard
else
    disp([mfilename,':: Derivatives of dynamic model wrt y(+0-), x and theta up to order ',sprintf('%0.0f',max_deriv_order),'. ',...
        sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])
end
%----------------------
routines.symbolic.probs_times_dynamic={original_funcs,wrt};

% % profile off, profile viewer
% % keyboard

% static model wrt y
%--------------------
static_incidence=zeros(orig_endo_nbr,3);
static_incidence(:,2)=1:orig_endo_nbr;
wrt=dynamic_differentiation_list(static_incidence,0,[]);

routines.static=utils.code.code2func(static.shadow_model);
[routines.static_derivatives,numEqtns,numVars,jac_toc,original_funcs]=...
    differentiate_system(...
    routines.static,dictionary.input_list,wrt,1);
routines.symbolic.static={original_funcs,wrt};
disp([mfilename,':: 1st-order derivatives of static model wrt y(0). ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])

% balanced growth path model wrt y
%---------------------------------
bgp_incidence=zeros(2*orig_endo_nbr,3);
bgp_incidence(:,2)=1:2*orig_endo_nbr;
wrt=dynamic_differentiation_list(bgp_incidence,0,[]);

routines.static_bgp=utils.code.code2func(static.shadow_BGP_model);
[routines.static_bgp_derivatives,numEqtns,numVars,jac_toc,original_funcs]=...
    differentiate_system(...
    routines.static_bgp,...
    dictionary.input_list,wrt,1);
routines.symbolic.static_bgp={original_funcs,wrt};
disp([mfilename,':: 1st-order derivatives of static BGP model wrt y(0). ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])
% dynamic model wrt param
%------------------------
param_nbr = numel(dictionary.parameters);
wrt=dynamic_differentiation_list([],0,[],1:param_nbr);
ppdd=@(x)x;%dynamic.shadow_model;
if ~dictionary.definitions_inserted
    if DefaultOptions.definitions_in_param_differentiation
        ppdd=@(x)parser.replace_definitions(x,shadow_definitions);
    else
        disp([mfilename,':: definitions not taken into account in the computation of derivatives wrt parameters'])
    end
end
[routines.parameter_derivatives,numEqtns,numVars,jac_toc,original_funcs]=...
    differentiate_system(...
    ppdd(dynamic.shadow_model),...
    dictionary.input_list,wrt,1);
routines.symbolic.parameters={original_funcs,wrt};
disp([mfilename,':: first-order derivatives of dynamic model wrt param. ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])

%% planner objective

if is_model_with_planner_objective
    % we use the same order as in the dynamic model to avoid a mismatch
    %------------------------------------------------------------------
    optimal_policy_incidence=dictionary.order_var.before_solve;
    wrt=dynamic_differentiation_list(optimal_policy_incidence,0,[]);
    
    routines.planner_objective=utils.code.code2func(dictionary.planner_system.shadow_model(1));
    [routines.planner_objective_derivatives,numEqtns,numVars,jac_toc,...
        original_funcs]=differentiate_system(...
        dictionary.planner_system.shadow_model(1),...
        dictionary.input_list,wrt,2);
    disp([mfilename,':: 1st and 2nd-order derivatives of planner objective wrt y(0). ',...
        sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.0f',jac_toc),' seconds'])
    % add the loss, the commitment degree and discount
    %-------------------------------------------------
    planner_shadow_model=strrep(strrep(dictionary.planner_system.shadow_model,'discount-',''),'commitment-','');
    routines.planner_loss_commitment_discount=utils.code.code2func(planner_shadow_model,dictionary.input_list);
    routines.symbolic.planner_objective={original_funcs,wrt};
end
%% Add final variables list to the dictionary
% the unsorted variables are variables sorted according to their order in
% during the solving of the model.
unsorted_endogenous=dictionary.orig_endogenous(dictionary.order_var.before_solve);
logical_incidence=dictionary.lead_lag_incidence.before_solve(dictionary.order_var.before_solve,:);

dictionary.NumberOfEquations=sum(equation_type==1);

% finally check that the number of equations is consistent with the number
% of variables

dictionary.is_sticky_information_model=ismember('sticky_information_lambda',{dictionary.parameters.name});

dictionary.is_hybrid_expectations_model=ismember('hybrid_expectations_lambda',{dictionary.parameters.name});

dictionary.is_dsge_var_model=ismember('dsge_prior_weight',{dictionary.parameters.name});

dictionary.is_optimal_policy_model=is_model_with_planner_objective && ...
    dictionary.NumberOfEquations<numel(dictionary.orig_endogenous);

dictionary.is_optimal_simple_rule_model=is_model_with_planner_objective && ...
    dictionary.NumberOfEquations==numel(dictionary.orig_endogenous);

if dictionary.is_sticky_information_model && dictionary.is_hybrid_expectations_model
    error([mfilename,':: you are not allowed to solve a model with both hybrid expectations and sticky information'])
end

dictionary.forward_looking_ids='NA';
dictionary.equations_reordering_for_multipliers=[];
unsorted_extra=unsorted_endogenous(1:0);
if dictionary.is_optimal_policy_model
    if dictionary.is_sticky_information_model
        error([mfilename,':: you are not allowed to solve a sticky information model with loose commitment'])
    end
    if dictionary.is_hybrid_expectations_model
        error([mfilename,':: you are not allowed to solve a hybrid expectations model with loose commitment'])
    end
    
    % lagrange multipliers
    %---------------------
    for eq=1:dictionary.NumberOfEquations
        new_var=struct('name',['mult_',sprintf('%0.0f',eq)],'tex_name','',...
            'max_lead',-equations_maxLag_maxLead(eq,1),... % lags govern the leads
            'max_lag',-equations_maxLag_maxLead(eq,2),... % leads govern the lags
            'is_log_var',false,'is_auxiliary',false,'is_trans_prob',false);
        for ifield=1:numel(useless_fields)
            new_var.(useless_fields{ifield})=nan;
        end
        unsorted_extra=[unsorted_extra,new_var];
    end
    LLI_EXTRA=[-equations_maxLag_maxLead(:,1),...
        ones(dictionary.NumberOfEquations,1),...
        -equations_maxLag_maxLead(:,2)];
    parts=utils.solve.partition_variables(LLI_EXTRA);
    unsorted_extra=unsorted_extra(parts.order_var);
    LLI_EXTRA=LLI_EXTRA(parts.order_var,:);
    dictionary.equations_reordering_for_multipliers=parts.order_var;
else
    added=0;
    assert(numel(dictionary.orig_endogenous)==sum(equation_type==1),...
        '# equations different from # endogenous variables')
    if dictionary.is_sticky_information_model
        dictionary.forward_looking_ids=find(dictionary.lead_lag_incidence.before_solve(:,3));
        for ii=1:numel(dictionary.forward_looking_ids)
            id=dictionary.forward_looking_ids(ii);
            new_var=struct('name',['SI_',dictionary.orig_endogenous{id}],'tex_name','',...
                'max_lead',0,'max_lag',0,'is_log_var',false,...
                'is_auxiliary',true,'is_trans_prob',false);
            for ifield=1:numel(useless_fields)
                new_var.(useless_fields{ifield})=nan;
            end
            unsorted_extra=[unsorted_extra,new_var];
        end
        added=numel(dictionary.forward_looking_ids);
    end
    % this will be used for the determination of the status of the variables.
    % the added variables are given status of static, which is misleading when
    % added>0
    LLI_EXTRA=repmat([0,1,0],added,1);
end
unsorted_endogenous=[unsorted_endogenous,unsorted_extra];
logical_incidence=[logical_incidence;LLI_EXTRA];

% now we can resort the final variables
[~,tags]=sort({unsorted_endogenous.name});
logical_incidence=logical_incidence(tags,:);
dictionary.endogenous=unsorted_endogenous(tags);

% update the lead-lag incidence and the order of the variables
%-------------------------------------------------------------
dictionary.lead_lag_incidence.after_solve=logical_incidence;
dictionary.lead_lag_incidence.after_solve(dictionary.lead_lag_incidence.after_solve~=0)=1:nnz(dictionary.lead_lag_incidence.after_solve);
% update the topology of the solution
%------------------------------------
[dictionary.siz.after_solve,...
    dictionary.locations.after_solve.t,...
    dictionary.locations.after_solve.z,~,...
    dictionary.order_var.after_solve,...
    dictionary.inv_order_var.after_solve]=...
    utils.solve.solution_topology(...
    dictionary.lead_lag_incidence.after_solve,...
    exo_nbr,... number of shocks
    0,... number of shocks periods beyond the current
    numel(switching_parameters_leads_index)... future switching parameters
    );
dictionary.switching_parameters_leads_index=switching_parameters_leads_index;

% reorder the solution of loose commitment or sticky information so that it
% conforms with order_var
%--------------------------------------------------------------------------
dictionary.reordering_index=locate_variables(...
    {dictionary.endogenous(dictionary.order_var.after_solve).name},{unsorted_endogenous.name});
dictionary.reordering_index=transpose(dictionary.reordering_index(:));
clear unsorted_endogenous

%% With all variables known, we can do the complementarity
[~,~,~,~,~,shadow_complementarity]=parser.shadowize(dictionary,mcp_eqtns,...
    mcp_type,orig_endogenous_current,overall_max_lead_lag);

routines.complementarity=utils.code.code2func(shadow_complementarity,parser.input_list);

%% give greek names to endogenous, exogenous, parameters
dictionary.endogenous=parser.greekify(dictionary.endogenous);
dictionary.orig_endogenous=parser.greekify(dictionary.orig_endogenous);
dictionary.exogenous=parser.greekify(dictionary.exogenous);
dictionary.parameters=parser.greekify(dictionary.parameters);

dictionary.is_purely_forward_looking_model=false;
dictionary.is_purely_backward_looking_model=false;
dictionary.is_hybrid_model=any(dictionary.lead_lag_incidence.before_solve(:,1)) && any(dictionary.lead_lag_incidence.before_solve(:,3));
if ~dictionary.is_hybrid_model
    if any(dictionary.lead_lag_incidence.before_solve(:,1))
        dictionary.is_purely_forward_looking_model=true;
    elseif any(dictionary.lead_lag_incidence.before_solve(:,3))
        dictionary.is_purely_backward_looking_model=true;
    end
end

% format endogenous, parameters, observables, etc
endogenous=dictionary.endogenous;
dictionary.endogenous=struct();
dictionary.endogenous.name={endogenous.name};
dictionary.endogenous.tex_name={endogenous.tex_name};
dictionary.endogenous.is_original=sparse(ismember(dictionary.endogenous.name,{dictionary.orig_endogenous.name}));
dictionary.endogenous.number=full([sum(dictionary.endogenous.is_original),numel(dictionary.endogenous.is_original)]);
dictionary.endogenous.is_lagrange_multiplier=sparse(strncmp('mult_',{endogenous.name},5));
dictionary.endogenous.is_static=sparse((~logical_incidence(:,1)&~logical_incidence(:,3))');
dictionary.endogenous.is_predetermined=sparse((~logical_incidence(:,1)&logical_incidence(:,3))');
dictionary.endogenous.is_pred_frwrd_looking=sparse((logical_incidence(:,1) & logical_incidence(:,3))');
dictionary.endogenous.is_state=dictionary.endogenous.is_predetermined|...
    dictionary.endogenous.is_pred_frwrd_looking;
dictionary.endogenous.is_frwrd_looking=sparse((logical_incidence(:,1) & ~logical_incidence(:,3))');
dictionary.endogenous.is_log_var=sparse([endogenous.is_log_var]);
dictionary.endogenous.is_log_expanded=sparse(false(size(dictionary.endogenous.is_log_var)));
dictionary.endogenous.is_auxiliary=sparse([endogenous.is_auxiliary]);
dictionary.endogenous.is_affect_trans_probs=sparse([endogenous.is_trans_prob]);
clear endogenous logical_incidence

exogenous=dictionary.exogenous;
dictionary.exogenous=struct();
dictionary.exogenous.name={exogenous.name};
dictionary.exogenous.tex_name={exogenous.tex_name};
dictionary.exogenous.is_observed=sparse(ismember(dictionary.exogenous.name,{dictionary.observables.name}));
dictionary.exogenous.number=full([sum(~dictionary.exogenous.is_observed),sum(dictionary.exogenous.is_observed)]);
dictionary.exogenous.is_in_use=sparse([exogenous.is_in_use]);
dictionary.exogenous.shock_horizon=sparse(1,sum(dictionary.exogenous.number));%double(~dictionary.exogenous.is_observed);
clear exogenous

observables=dictionary.observables;
dictionary.observables=struct();
dictionary.observables.name={observables.name};
dictionary.observables.is_endogenous=sparse(ismember(dictionary.observables.name,dictionary.endogenous.name));
tex_names={observables.tex_name};
state_endo=locate_variables(dictionary.observables.name,dictionary.endogenous.name,true);
state_endo(isnan(state_endo))=0;
tex_names(state_endo>0)=dictionary.endogenous.tex_name(nonzeros(state_endo));
state_exo=locate_variables(dictionary.observables.name,dictionary.exogenous.name,true);
state_exo(isnan(state_exo))=0;
tex_names(state_exo>0)=dictionary.exogenous.tex_name(nonzeros(state_exo));
state_id=state_endo+state_exo*1i;
dictionary.observables.state_id=state_id(:).';
dictionary.observables.tex_name=tex_names;
dictionary.observables.number=full([sum(dictionary.observables.is_endogenous),sum(~dictionary.observables.is_endogenous)]);
clear observables

parameters=dictionary.parameters;
dictionary.parameters=struct();
dictionary.parameters.name={parameters.name};
dictionary.parameters.tex_name={parameters.tex_name};
dictionary.parameters.is_switching=sparse([parameters.is_switching]);
dictionary.parameters.is_trans_prob=sparse([parameters.is_trans_prob]);
dictionary.parameters.is_measurement_error=sparse([parameters.is_measurement_error]);
% after taking some sparse above, we have to make sure that the matrix
% below is full
dictionary.parameters.number=full([sum(~dictionary.parameters.is_switching),sum(dictionary.parameters.is_switching)]);
dictionary.parameters.is_in_use=sparse([parameters.is_in_use]);
dictionary.parameters.governing_chain=[parameters.governing_chain];
clear parameters

% variable names for the original endogenous but without the augmentation
% add-ons
dictionary.orig_endo_names_current={orig_endogenous_current.name};

% equations
%----------
dictionary.equations.dynamic=dynamic.model;
dictionary.equations.shadow_dynamic=dynamic.shadow_model;
dictionary.equations.static=static.model;
dictionary.equations.shadow_static=static.shadow_model;
dictionary.equations.shadow_balanced_growth_path=static.shadow_BGP_model;
dictionary.equations.number=numel(dynamic.model);
dictionary.is_imposed_steady_state=static.is_imposed_steady_state;
dictionary.is_unique_steady_state=static.is_unique_steady_state;

% dictionary.model_derivatives=dynamic.model_derivatives;
% dictionary.steady_state_shadow_model=static.steady_state_shadow_model;
%  char({dictionary.equations.shadow_balanced_growth_path})
definitions=dictionary.definitions;
dictionary.definitions=struct();
dictionary.definitions.dynamic={definitions.model};
dictionary.definitions.dynamic=dictionary.definitions.dynamic(:);
dictionary.definitions.shadow_dynamic={definitions.shadow};
dictionary.definitions.shadow_dynamic=dictionary.definitions.shadow_dynamic(:);
dictionary.definitions.number=numel(definitions);
clear definitions static dynamic
dictionary.routines=routines;

dictionary=orderfields(dictionary);

end

function [derivs,numEqtns,numVars,jac_toc,original_funcs]=differentiate_system(myfunc,input_list,wrt,order)

numEqtns=numel(myfunc);

myfunc=parser.remove_handles(myfunc);

myfunc=parser.symbolic_model(myfunc,input_list);

% list of symbols
symb_list=parser.collect_symbolic_list(myfunc,strcat(input_list,'_'));
% force 's0' and 's1' to enter the list
state_inputs={'s0','s1'};
input_list=input_list(:)';
for ii=1:numel(state_inputs)
    if ~any(strcmp(symb_list,state_inputs{ii}))
        symb_list=[symb_list,state_inputs{ii}];
    end
end
% sorting will be useful if we need to debug
symb_list=sort(symb_list);

args=splanar.initialize(symb_list,wrt);
numVars=numel(wrt);
original_funcs=myfunc;
for ifunc=1:numEqtns
    [occur,myfunc{ifunc}]=parser.find_occurrences(myfunc{ifunc},symb_list);
    % re-create the function
    var_occur=symb_list(occur);
    argfun=cell2mat(strcat(var_occur,','));
    myfunc{ifunc}=str2func(['@(',argfun(1:end-1),')',myfunc{ifunc}]);
    original_funcs{ifunc}=myfunc{ifunc};
    arg_occur=args(occur);
    myfunc{ifunc}=myfunc{ifunc}(arg_occur{:});
end

verbose=false;
tic
derivs=splanar.differentiate(myfunc,numVars,order,verbose);
derivs=splanar.print(derivs,input_list);
jac_toc=toc;

end

function [wrt,v,locations,siz,order_var,inv_order_var,steady_state_index]=...
    dynamic_differentiation_list(LLI,exo_nbr,spindex,pindex)%% partition the endogenous
if nargin<4
    pindex=[];
end

% order the list of the variables to differentiate according to
%---------------------------------------------------------------
v={
    'b_plus',0
    'f_plus',0
    's_0',0
    'p_0',0
    'b_0',0
    'f_0',0
    'p_minus',0
    'b_minus',0
    'e_0',0
    'theta_plus',0
    };
fields=v(:,1);

locations=struct();

[siz,locations.t,locations.z,~,order_var,inv_order_var]=...
    utils.solve.solution_topology(...
    LLI,...
    exo_nbr,... number of shocks
    0,... number of shocks periods beyond the current
    numel(spindex)... future switching parameters
    );
v{strcmp(v(:,1),'s_0'),2}=siz.ns;

v(strncmp(v(:,1),'p',1),2)={siz.np};

v(strncmp(v(:,1),'b',1),2)={siz.nb};

v(strncmp(v(:,1),'f',1),2)={siz.nf};

steady_state_index=struct();
if isempty(LLI)
    ywrt={};
else
    % endogenous
    %-----------
    llio=LLI(order_var,:);
    yindex=nonzeros(llio(:))';
    ywrt=cellstr(strcat('y_',num2str(yindex(:))));
    steady_state_index.y=order_var([find(llio(:,1))',find(llio(:,2))',find(llio(:,3))']); %<-- yindex=[LLI([both,frwrd],1)',LLI(order_var,2)',LLI([pred,both],3)];
    ywrt=cellfun(@(x)x(~isspace(x)),ywrt,'uniformOutput',false);
end
% exogenous
%----------
if exo_nbr
    xindex=1:exo_nbr;
    xwrt=cellstr(strcat('x_',num2str(xindex(:))));
    xwrt=cellfun(@(x)x(~isspace(x)),xwrt,'uniformOutput',false);
    steady_state_index.x=xindex(:)';
else
    xwrt={};
end
v{strcmp(v(:,1),'e_0'),2}=siz.ne;

% switching parameters
%---------------------
spwrt=process_parameters(spindex,'s');
steady_state_index.theta_plus=spindex(:)';
v{strcmp(v(:,1),'theta_plus'),2}=siz.ntheta;

% constant parameters
%--------------------
pwrt=process_parameters(pindex,'');

% differentiation list
%---------------------
wrt=[ywrt(:)',xwrt(:)',spwrt(:)',pwrt(:)'];
steady_state_index.wrt=wrt;

% v-locations
%------------
numbers=cell2mat(v(:,2));
for ifield=1:numel(fields)
    ff=fields{ifield};
    underscore=strfind(ff,'_');
    ff2=ff(1:underscore(1)-1);
    locations.v.(ff)=sum(numbers(1:ifield-1))+(1:siz.(['n',ff2]));
end
siz.nv=sum(numbers);
locations.v.bf_plus=[locations.v.b_plus,locations.v.f_plus];
locations.v.pb_minus=[locations.v.p_minus,locations.v.b_minus];
locations.v.t_0=siz.nb+siz.nf+(1:siz.ns+siz.np+siz.nb+siz.nf); % =[locations.v.s_0,locations.v.p_0,locations.v.b_0,locations.v.f_0];

    function pwrt=process_parameters(index,prefix)
        if isempty(index)
            pwrt={};
        else
            if ischar(index)
                index=cellstr(index);
            end
            pwrt=cellstr(strcat(prefix,'param_',num2str(index(:))));
            pwrt=cellfun(@(x)x(~isspace(x)),pwrt,'uniformOutput',false);
        end
    end
end
