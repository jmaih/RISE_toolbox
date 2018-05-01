function dictionary=parse(FileName,varargin)
% parse -- parser for dsge models
%
% ::
%
%
%   dictionary=parse(FileName,varargin)
%
% Args:
%
%    - **FileName** [char|cellstr]: if "char", name of the model file. The
%    file should have extensions rs, rz or dsge. If "cellstr" each cell
%    contains the name of a separate model file. The files are then meant to
%    be combined into one single model.
%
%    - **varargin** []: pairwise arguments with possiblities as follows:
%      - **parameter_differentiation** [true|{false}]: compute or not
%      parameter derivatives
%      - **definitions_inserted** [true|{false}]: substitute definitions
%      - **definitions_in_param_differentiation** [true|{false}]: insert or
%      not definitions in equations before differentiating with respect to
%      parameters
%      - **saveas** [true|false|char|{''}]: save the possibly expanded model
%      file. If "true", the name of the main original file is used appended with
%      "_expanded.dsge". Alternatively, the user can provide a name under
%      which he wants the file to be saved.
%      - **max_deriv_order** [integer|{1}]: order for symbolic
%      differentiation. It is recommended to set to 1, especially for large
%      models in case one does not intend to solve higher-order approximations
%      - **parse_debug** [true|{false}]: debugging in the parser
%      - **add_welfare** [true|{false}]: add the welfare equation when doing
%      optimal policy. N.B: The welfare variable, WELF is the true welfare
%      multiplied by (1-discount). The within-period utility variable, UTIL is
%      automatically added. The reason why welfare is not automatically added
%      is that oftentimes models with that equation do not solve.
%      - **rise_flags** [struct|cell]: instructions for the partial parsing of
%      the rise file. In case of a cell, the cell should be a k x 2 cell,
%      where the first column collects the conditional parsing names and the
%      second column the values.
%
% Returns:
%    :
%
%    - **dictionary** [struct]: elements need by the dsge class for
%    constructing an instance.
%
% Note:
%
%    - In RISE it is possible to declare exogenous and make them observable at
%    the same time. The exogenous that are observed are determisitic. This
%    is the way to introduce e.g. time trends. This strategy also opens the
%    door for estimating partial equilibrium models
%
% Example:
%
%    See also:

DefaultOptions=...
    struct('definitions_in_param_differentiation',false,...
    'saveas','',...
    'max_deriv_order',1,'definitions_inserted',false,...
    'parse_debug',false,...
    'add_welfare',false,...
    'parameter_differentiation',false);
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
dictionary.add_welfare=DefaultOptions.add_welfare;
parameter_differentiation=DefaultOptions.parameter_differentiation;

%% set various blocks

% first output: the dictionary.filename
valid_extensions={'.rs','.rz','.dsge'};
VE=parser.cell2matize(valid_extensions);
FileName=regexprep(FileName,'\s*','');
FileName=regexp(FileName,['(?<fname>\w+[^\.]*)(?<ext>',VE,'?)'],'names');
if iscell(FileName)
    FileName=[FileName{:}];
end
for id=1:numel(FileName)
    if isempty(FileName(id).ext)
        iext=0; found=false;
        while ~found && iext<numel(valid_extensions)
            iext=iext+1;
            found=exist([FileName(id).fname,valid_extensions{iext}],'file');
        end
        if ~found
            error([mfilename,':: ',FileName(id).fname,'.rs or ',...
                FileName(id).fname,'.rz  or ',FileName(id).fname,'.dsge not found'])
        end
        FileName(id).ext=valid_extensions{iext};
    else
        if ~exist([FileName(id).fname,FileName(id).ext],'file')
            error([mfilename,':: ',[FileName(id).fname,FileName(id).ext],' not found'])
        end
    end
end

filename_=strrep(parser.cell2matize({FileName.fname}),'|','+');
filename_=strrep(filename_,'(','');
filename_=strrep(filename_,')','');

dictionary.filename=filename_;

% read file and remove comments
% RawFile=read_file(FileName,DefaultOptions.rise_flags);

RawFile=parser.preparse(FileName,DefaultOptions.rise_flags);

logic=islogical(DefaultOptions.saveas) && DefaultOptions.saveas;

hasname= ~isempty(DefaultOptions.saveas) && ischar(DefaultOptions.saveas);

newname='';
if logic
    newname=[strrep(filename_,'+','_'),'_expanded.dsge'];
elseif hasname
    newname=DefaultOptions.saveas;
    thedot=strfind(newname,'.');
    if isempty(thedot) %#ok<STREMP>
        newname=[newname,'.dsge'];
    end
end
% write the expanded version
if ~isempty(newname)
    newfile='';
    fid=fopen(newname,'w');
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

old_endo_names={dictionary.endogenous.name};
logvarnames={dictionary.log_vars.name};
for ivar=1:numel(logvarnames)
    loc=strcmp(logvarnames{ivar},{dictionary.endogenous.name});
    dictionary.endogenous(loc).is_log_var=true;
end
%% Model block
% now with the endogenous, exogenous, parameters in hand, we can process

[Model_block,dictionary,blocks]=parser.parse_model(dictionary,blocks);

%% optimal policy and optimal simple rule block

[dictionary,blocks,Model_block,PlannerObjective_block,jac_toc_]=...
            parser.parse_optimal_policy(Model_block,dictionary,blocks);        
        
if ~isempty(jac_toc_)
    disp([mfilename,':: First-order conditions of optimal policy :',sprintf('%0.4f',jac_toc_),' seconds'])
end

%% after parsing the model block, update the markov chains (time-varying probabilities)
% sorting the endogenous switching probabilities is more or less useless
dictionary.time_varying_probabilities=sort(dictionary.time_varying_probabilities);

%% Now we can re-order the original (and augmented) endogenous
[~,tags]=sort({dictionary.endogenous.name});
dictionary.endogenous=dictionary.endogenous(tags);
% also renew the endogenous quick list for status determination under
% shadowization, where the precise locations now matter!!!
%----------------------------------------------------------------------
dictionary.endogenous_list={dictionary.endogenous.name};
% no need to renew the parameters and the exogenous since new ones have not
% been created

%% Now we re-write the model and update leads and lags
% at the same time, construct the incidence and occurrence matrices
orig_endo_nbr=numel(dictionary.endogenous);

[equation_type,occurrence,fast_sstate_occurrence]=equation_types_and_variables_occurrences();

%% Steady state Model block
% now with the endogenous, exogenous, parameters in hand, we can process
% the steady state model block
[static,SteadyStateModel_block,auxiliary_steady_state_equations,...
    dictionary,blocks]=parser.parse_steady_state(dictionary,blocks);
%% exogenous definitions block

[dictionary,blocks]=parser.parse_exogenous_definitions(dictionary,blocks);

%% parameterization block

do_parameterization();

%% parameter restrictions block.

do_parameter_restrictions();

%% Lump together the model and steady-state model
static_mult_equations=[];
utils_derivatives=[];
if dictionary.is_optimal_policy_model
    
    static_mult_equations=dictionary.planner_system.static_mult_equations{1};
    
end

if dictionary.is_model_with_planner_objective

    utils_derivatives=dictionary.planner_system.utils_derivs.derivatives;
    
end

AllModels=[Model_block
    SteadyStateModel_block
    auxiliary_steady_state_equations
    PlannerObjective_block
    utils_derivatives
    static_mult_equations];
% steady state equations (including auxiliary) are identified by number 5
aux_ss_eq_nbr=size(auxiliary_steady_state_equations,1);
ss_eq_nbr=size(SteadyStateModel_block,1);
% dictionary.planner_system objective equations are identified by number 6
planobj_eq_nbr=size(PlannerObjective_block,1);
osr_derivs_eqn_nbr=size(utils_derivatives,1);
% mult steady state identified by number 7
stat_mult_eq_nbr=size(static_mult_equations,1);
equation_type=[equation_type
    5*ones(ss_eq_nbr+aux_ss_eq_nbr,1)
    6*ones(planobj_eq_nbr+osr_derivs_eqn_nbr,1)
    7*ones(stat_mult_eq_nbr,1)];

clear Model_block SteadyStateModel_block PlannerObjective_block
%% Incidence matrix
% Now and only now can we build the incidence matrix
do_incidence_matrix();

%% models in shadow/technical/tactical form

[dictionary,...
    dynamic,...
    stat,...
    defs,...
    shadow_tvp,...
    shadow_complementarity]=parser.shadowize(dictionary,AllModels,...
    equation_type);

% occurrence of parameters in case of swap when solving steady state
%--------------------------------------------------------------------
fast_sstate_occurrence_params=parser.occurrence_map(stat.shadow_fast_ssmodel,...
    'param',numel(dictionary.parameters),true);

% Remove the fast steady state from the dynamic!!!
% then build the bgp/static system
%---------------------------------------------------

orig_definitions=defs.original;
shadow_definitions=defs.shadow;

%% steady state and static models
static=utils.miscellaneous.mergestructures(static,stat);

static.shadow_steady_state_model=parser.substitute_definitions(...
    static.shadow_steady_state_model,shadow_definitions);

if dictionary.is_optimal_policy_model
    dictionary.planner_system.static_mult_equations{1}=...
        static.shadow_steady_state_model(ss_eq_nbr+aux_ss_eq_nbr+1:end);
end

static.shadow_steady_state_auxiliary_eqtns=...
    static.shadow_steady_state_model(ss_eq_nbr+(1:aux_ss_eq_nbr));

static.shadow_steady_state_model(ss_eq_nbr+1:end)=[];

dictionary.is_param_changed_in_ssmodel=parser.parameters_changed_in_ssmodel(...
    static.shadow_steady_state_model,'param',numel(dictionary.parameters));

static.shadow_steady_state_model=strrep(static.shadow_steady_state_model,...
    'x1_=','[x1_,fval,exitflag]=');
static.shadow_steady_state_model=strrep(static.shadow_steady_state_model,...
    ',x0_,',',x0_,options,');
old_shadow_steady_state_model=static.shadow_steady_state_model;
static.shadow_steady_state_model=cell(0,1);
fsolve_nbr=0;
for ii=1:numel(old_shadow_steady_state_model)
    eq_i=old_shadow_steady_state_model{ii};
    if ~isempty(strfind(eq_i,'argzero')) %#ok<STREMP>
        eq_i=strrep(eq_i,'argzero','fsolve');
        static.shadow_steady_state_model=[static.shadow_steady_state_model;{eq_i}];
        eq_i={'retcode=1-(exitflag==1);'};
        if ii<numel(old_shadow_steady_state_model)
            eq_i=[eq_i;{'if ~retcode,'}]; %#ok<*AGROW>
            fsolve_nbr=fsolve_nbr+1;
        end
        static.shadow_steady_state_model=[static.shadow_steady_state_model;eq_i];
    else
        static.shadow_steady_state_model=[static.shadow_steady_state_model;{eq_i}];
    end
end
for ii=1:fsolve_nbr
    static.shadow_steady_state_model=[static.shadow_steady_state_model;{'end;'}];
end
clear old_shadow_steady_state_model
% now add the prelude
if ~isempty(static.shadow_steady_state_model)
    static.shadow_steady_state_model=[{['if ~exist(''y'',''var''),y=zeros(',sprintf('%0.0f',orig_endo_nbr),',1); end;']};static.shadow_steady_state_model];
end
%% replace the list of definition names with definition equations
dictionary.definitions=struct('name',dictionary.definitions(:),...
    'model',orig_definitions(:),'shadow',shadow_definitions(:));

%% Routines structure
routines=struct();

% With all variables known, we can do the complementarity
%---------------------------------------------------------
routines.complementarity=utils.code.code2func(shadow_complementarity,parser.input_list);

% definitions routine
%---------------------
routines.definitions=utils.code.code2func(...
    parser.substitute_definitions(shadow_definitions),'param');

% steady state model routine (cannot be written as a function)
%--------------------------------------------------------------
routines.steady_state_model=parser.substitute_definitions(...
    static.shadow_steady_state_model,shadow_definitions);
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

% v=[f+,b+,s0,p_0,b_0,f_0,p_minus,b_minus,e_0]
routines.dynamic=utils.code.code2func(dynamic.shadow_model);

[wrt,dictionary.v,...
    dictionary.locations.before_solve,...
    dictionary.siz.before_solve,...
    dictionary.order_var,...
    dictionary.inv_order_var,...
    dictionary.steady_state_index]=dynamic_differentiation_list(...
    dictionary.lead_lag_incidence.before_solve,exo_nbr,[]);
%----------------------
if dictionary.parse_debug
    
    profile off
    
    profile on
    
end

routines.probs_times_dynamic=parser.burry_probabilities(dynamic.shadow_model,myifelseif);

[routines.probs_times_dynamic_derivatives,numEqtns,numVars,jac_toc,...
    original_funcs]=parser.differentiate_system(...
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

% The so-called static model may also have some growth in which case the
% relevant thing to differentiate is not the static model...
% static model wrt y
%--------------------
% static_incidence=zeros(orig_endo_nbr,3);
% static_incidence(:,2)=1:orig_endo_nbr;
% wrt=dynamic_differentiation_list(static_incidence,0);

% [routines.static_derivatives,numEqtns,numVars,jac_toc,original_funcs]=...
%     differentiate_system(routines.static,dictionary.input_list,wrt,1);
% routines.symbolic.static={original_funcs,wrt};
% disp([mfilename,':: 1st-order derivatives of static model wrt y(0). ',...
%     sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])
% 
% dynamic model wrt param
%------------------------
if parameter_differentiation
    
    param_nbr = numel(dictionary.parameters);
    
    wrt=dynamic_differentiation_list([],0,1:param_nbr);
    
    ppdd=@(x)x;%dynamic.shadow_model;
    
    if ~dictionary.definitions_inserted
    
        if DefaultOptions.definitions_in_param_differentiation
        
            ppdd=@(x)parser.replace_definitions(x,shadow_definitions);
        
        else
            
            disp([mfilename,':: definitions not taken into account in the computation of derivatives wrt parameters'])
        
        end
        
    end
    
    [routines.parameter_derivatives,numEqtns,numVars,jac_toc,original_funcs]=...
        parser.differentiate_system(...
        ppdd(dynamic.shadow_model),...
        dictionary.input_list,wrt,1);
    
    routines.symbolic.parameters={original_funcs,wrt};
    
    disp([mfilename,':: first-order derivatives of dynamic model wrt param. ',...
        sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])

end
%% optimal policy and optimal simple rules routines
%-----------------------------------------
if dictionary.is_model_with_planner_objective
    planner_shadow_model=strrep(strrep(dictionary.planner_system.shadow_model,'discount-',''),'commitment-','');
    
    if dictionary.is_optimal_policy_model
        tmp=parser.replace_steady_state_call(dictionary.planner_system.static_mult_equations{1});
        tmp=utils.code.code2func(tmp,parser.input_list);
        
        routines.planner_static_mult=tmp;
        routines.planner_static_mult_support=...
            dictionary.planner_system.static_mult_equations(2:end);
        dictionary.planner_system=rmfield(dictionary.planner_system,...
            'static_mult_equations');
    end
    osr_=dictionary.planner_system.utils_derivs;
    endo_names={dictionary.endogenous.name};
    ordered_endo_names=endo_names(dictionary.order_var);
    der_reo=locate_variables(osr_.wrt,ordered_endo_names);
    routines.planner_osr_support=struct('derivatives_re_order',der_reo,...
        'partitions',osr_.partitions,'nwrt',osr_.nwrt,...
        'map',vec(cell2mat(osr_.map(:,2).')),'size',osr_.size);
    % we take the second column since the first column with the
    % equation numbers do not matter: originally we had only one equation
    clear osr_
    % add the loss, the commitment degree and discount
    %-------------------------------------------------
    routines.planner_loss_commitment_discount=utils.code.code2func(planner_shadow_model,dictionary.input_list);
    
    routines.planner_objective=utils.code.code2func(dictionary.planner_system.shadow_model(1));

    dictionary.planner_system=rmfield(dictionary.planner_system,...
        {'shadow_model','model'});
end
%% Add final variables list to the dictionary
% the unsorted variables are variables sorted according to their order in
% during the solving of the model.
unsorted_endogenous=dictionary.endogenous(dictionary.order_var);
logical_incidence=dictionary.lead_lag_incidence.before_solve(dictionary.order_var,:);

dictionary.NumberOfEquations=sum(equation_type==1);

% now we can resort the final variables
[~,tags]=sort({unsorted_endogenous.name});
logical_incidence=logical_incidence(tags,:);
dictionary.endogenous=unsorted_endogenous(tags);

% update the lead-lag incidence and the order of the variables: with the
% new settings, this is not expected to ever change
%--------------------------------------------------------------------------
dictionary.lead_lag_incidence.after_solve=logical_incidence;
dictionary.lead_lag_incidence.after_solve(dictionary.lead_lag_incidence.after_solve~=0)=1:nnz(dictionary.lead_lag_incidence.after_solve);

% update the topology of the solution: with the new settings, this is not
% expected to ever change
%--------------------------------------------------------------------------
[dictionary.siz.after_solve,...
    dictionary.locations.after_solve.t,...
    dictionary.locations.after_solve.z]=...
    utils.solve.solution_topology(...
    dictionary.lead_lag_incidence.after_solve,...
    exo_nbr,... number of shocks
    0); % number of shocks periods beyond the current

clear unsorted_endogenous

%% clean up
dictionary.routines=routines;

dictionary.occurrence=occurrence;

dictionary.fast_sstate_occurrence=struct('y',fast_sstate_occurrence,...
    'param',fast_sstate_occurrence_params);

dictionary=parser.dictionary_cleanup(dictionary,dynamic,static,old_endo_names,logical_incidence);

    function do_incidence_matrix()
        
        before_solve=zeros(3,orig_endo_nbr);
        for iii=1:3
            for jj=1:orig_endo_nbr
                if any(occurrence(:,jj,iii))
                    before_solve(iii,jj)=1;
                end
            end
        end
        before_solve=transpose(flipud(before_solve));
        before_solve(before_solve>0)=1:nnz(before_solve);
        
        appear_as_current=before_solve(:,2)>0;
        if any(~appear_as_current)
            disp('The following variables::')
            allendo={dictionary.endogenous.name};
            disp(allendo(~appear_as_current))
            error('do not appear as current')
        end
        dictionary.lead_lag_incidence.before_solve=before_solve;
    end

    function [equation_type,occurrence,fast_sstate_occurrence]=equation_types_and_variables_occurrences()
        number_of_equations=size(Model_block,1);
        occurrence=false(number_of_equations,orig_endo_nbr,3);
        fast_sstate_occurrence=occurrence;
        has_fast_sstate=cellfun(@(x)~isempty(x),Model_block(:,end));
        equation_type=ones(number_of_equations,1);
        for iii=1:number_of_equations
            % main equation
            eq_i_= Model_block{iii,1};
            run_occurrence(eq_i_);
            % bgp/sstate equation if not empty
            eq_i_= Model_block{iii,end};
            run_occurrence(eq_i_,false);
        end
        % replace the locations without steady state equations with the
        % normal equations
        %---------------------------------------------------------------
        fast_sstate_occurrence(~has_fast_sstate,:,:)=occurrence(~has_fast_sstate,:,:);
        % keep only the structural equations
        %------------------------------------
        occurrence=occurrence(equation_type==1,:,:);
        fast_sstate_occurrence=fast_sstate_occurrence(equation_type==1,:,:);
        neqtns=size(occurrence,1);
        nvars=size(occurrence,2);
        if neqtns>nvars
            error(['more equations (',int2str(neqtns),') than variables(',int2str(nvars),')'])
        elseif neqtns<nvars
            error(['more variables (',int2str(nvars),') than equations(',int2str(neqtns),')'])
        end
        
        function run_occurrence(eq_i_,is_dynamic)
            if nargin<2
                is_dynamic=true;
            end
            if isempty(eq_i_)
                return
            end
            if strcmp(Model_block{iii,4},'def') % <---ismember(eq_i_{1,1},dictionary.definitions) && strcmp(eq_i_{1,2}(1),'=')
                equation_type(iii)=2;
            elseif strcmp(Model_block{iii,4},'tvp') % <---ismember(eq_i_{1,1},dictionary.time_varying_probabilities) && strcmp(eq_i_{1,2}(1),'=')
                equation_type(iii)=3;
            elseif strcmp(Model_block{iii,4},'mcp') % <---ismember(eq_i_{1,1},dictionary.time_varying_probabilities) && strcmp(eq_i_{1,2}(1),'=')
                equation_type(iii)=4;
            end
            for ii2=1:size(eq_i_,2)
                if ~isempty(eq_i_{2,ii2})
                    var_loc=strcmp(eq_i_{1,ii2},{dictionary.endogenous.name});
                    lag_or_lead=eq_i_{2,ii2}+2;
                    if any(var_loc) && equation_type(iii)==2
                        error([mfilename,':: equation (',sprintf('%0.0f',iii),') detected to be a definition cannot contain variables'])
                    elseif equation_type(iii)==3 && ismember(lag_or_lead,[3,1])
                        error([mfilename,':: equation (',sprintf('%0.0f',iii),') detected to describe endogenous switching cannot contain leads or lags'])
                    end
                    if is_dynamic
                        occurrence(iii,var_loc,lag_or_lead)=true;
                    else
                        fast_sstate_occurrence(iii,var_loc,lag_or_lead)=true;
                    end
                end
            end
        end
    end

    function do_parameter_restrictions()
        
        current_block_id=find(strcmp('parameter_restrictions',{blocks.name}));
        
        [Param_rest_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'parameter_restrictions');
        
        dictionary.Param_rest_block=struct('original',...
            {blocks(current_block_id).listing(:,2)},...
            'parsed',{Param_rest_block(:,1)});
        % remove item from block
        blocks(current_block_id)=[];
    end

    function do_parameterization()
        
        current_block_id=find(strcmp('parameterization',{blocks.name}));
        
        dictionary.Parameterization_block=parser.capture_parameterization(dictionary,blocks(current_block_id).listing);
        
        % remove item from block
        blocks(current_block_id)=[];
    end
end


function [wrt,v,locations,siz,order_var,inv_order_var,steady_state_index]=...
    dynamic_differentiation_list(LLI,exo_nbr,pindex)%% partition the endogenous
if nargin<3
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
    };
fields=v(:,1);

locations=struct();

[siz,locations.t,locations.z,~,order_var,inv_order_var]=...
    utils.solve.solution_topology(...
    LLI,...
    exo_nbr,... number of shocks
    0);% number of shocks periods beyond the current
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

% constant parameters
%--------------------
pwrt=process_parameters(pindex,'');

% differentiation list
%---------------------
wrt=[ywrt(:)',xwrt(:)',pwrt(:)'];
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
