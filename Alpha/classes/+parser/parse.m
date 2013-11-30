function dictionary=parse(FileName,varargin)
% by the way, I can still declare exogenous and make them observable at the
% same time. The exogenous that are observed are determisitic. This opens
% the door for estimating partial equilibrium models

DefaultOptions=...
    struct('definitions_in_param_differentiation',true,...
    'rise_flags',struct(),'rise_save_macro',false);
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
    else
        error([mfilename,':: ',FileName,'.rs or ',FileName,'.rz not found'])
    end
else
    ext=FileName(loc:end);
%     if ~(strcmp(ext,'.rs')||strcmp(ext,'.dyn')||strcmp(ext,'.mod'))
    if ~ismember(ext,{'.rs','.rz'})%||strcmp(ext,'.dyn')||strcmp(ext,'.mod'))
        error([mfilename,':: Input file is expected to be a .rs or a .rz file'])
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
    fid=fopen([FileName(1:thedot-1),'_expanded.rs'],'w');
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
[Model_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'model');

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
auxiliary_steady_state_equations=cell(0,3);
orig_endogenous_current=dictionary.orig_endogenous; % these are the variables without the augmentation add-on
all_fields=fieldnames(parser.listing);
main_endo_fields={'name','tex_name','max_lead','max_lag','is_log_var','is_auxiliary'};
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
        new_var=struct('name',[vname,'_AUX_F_',sprintf('%0.0f',i2-1)],...
            'tex_name','','max_lead',1,'max_lag',0,'is_log_var',false,'is_auxiliary',true);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vold,1}',{';',[]}'],0,1}]; %#ok<*AGROW> %
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{vold,0}',{';',[]}'],0,0}]; %<-- add maxLag and maxLead
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
            'max_lag',-1,'is_log_var',false,'is_auxiliary',true);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vname,0}',{';',[]}'],0,0}]; %
        % The steady state is computed with zero shocks. and so instead of
        % setting the name of the shock, I write 0
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{'0',0}',{';',[]}'],0,0}]; %<-- add maxLag and maxLead
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
        new_var=struct('name',[vname,'_AUX_L_',sprintf('%0.0f',i2-1)],...
            'tex_name','','max_lead',0,'max_lag',-1,'is_log_var',false,...
            'is_auxiliary',true);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vold,-1}',{';',[]}'],-1,0}]; %
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{vold,0}',{';',[]}'],0,0}]; %<-- add maxLag and maxLead
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
    if ismember(eq_i{1,1},dictionary.definitions) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=2;
    elseif ismember(eq_i{1,1},dictionary.time_varying_probabilities) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=3;
    end
    for i2=1:size(eq_i,2)
        if ~isempty(eq_i{2,i2})
            vname=eq_i{1,i2};
            status=parser.determine_status(dictionary,vname);
            time=-1; new_item=false;
            if strcmp(status,'x')&& abs(eq_i{2,i2})>0
                new_item=true;
                if abs(eq_i{2,i2})==1
                    % change the name, but not the lag structure
                    eq_i{1,i2}=[vname,'_0'];
                else
                    % change both the name and the lag structure
                    eq_i{1,i2}=[vname,'_0','_AUX_L_',sprintf('%0.0f',abs(eq_i{2,i2})-1)];
                end
            elseif strcmp(status,'y')&& abs(eq_i{2,i2})>1
                new_item=true;
                % change both the name and the lag structure
                if eq_i{2,i2}>0
                    eq_i{1,i2}=[vname,'_AUX_F_',sprintf('%0.0f',abs(eq_i{2,i2})-1)];
                    time=1;
                elseif eq_i{2,i2}<0
                    eq_i{1,i2}=[vname,'_AUX_L_',sprintf('%0.0f',abs(eq_i{2,i2})-1)];
                end
            end
            if new_item
                eq_i{2,i2}=time;
            end
            var_loc=strcmp(eq_i{1,i2},{dictionary.orig_endogenous.name});
            lag_or_lead=eq_i{2,i2}+2;
            if equation_type(ii)==2
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
[SteadyStateModel_block,dictionary,static]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'steady_state_model',static);

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
[ExogenousDefinition_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'exogenous_definition');

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
% steady state equations are identified by number 4
planobj_eq_nbr=size(PlannerObjective_block,1);
% dictionary.planner_system objective equations are identified by number 5
equation_type=[equation_type;4*ones(ss_eq_nbr,1);5*ones(planobj_eq_nbr,1)];

% collect information about leads and lags which will be used to determine
% the status of the lagrange multipliers later on. But keep only the
% information about the structural equations
%--------------------------------------------------------------------------
equations_maxLag_maxLead=cell2mat(Model_block(equation_type==1,2:end));

clear Model_block SteadyStateModel_block PlannerObjective_block
%% Incidence matrix
% Now and only now can we build the incidence matrix
dictionary.Lead_lag_incidence=zeros(3,orig_endo_nbr);
for ii=1:3
    for jj=1:orig_endo_nbr
        if any(Occurrence(:,jj,ii))
            dictionary.Lead_lag_incidence(ii,jj)=1;
        end
    end
end
dictionary.Lead_lag_incidence=transpose(flipud(dictionary.Lead_lag_incidence));
dictionary.Lead_lag_incidence(dictionary.Lead_lag_incidence>0)=1:nnz(dictionary.Lead_lag_incidence);

appear_as_current=find(dictionary.Lead_lag_incidence(:,2));
if any(~appear_as_current)
    disp('The following variables::')
    disp(old_endo_names(~appear_as_current))
    error('do not appear as current')
end

%%

dynamic.model=cell(0,1);
dynamic.shadow_model=cell(0,1);
static.model=cell(0,1);
static.shadow_model=cell(0,1);
static.shadow_BGP_model=cell(0,1);
dictionary.planner_system.shadow_model=cell(0,1);

static.steady_state_shadow_model=cell(0,1);
shadow_definitions=cell(0,1);
orig_definitions=cell(0,1);
shadow_tvp=cell(0,1);

for ii=1:numel(equation_type)
    % we don't need the semicolons anymore
    eq_i=AllModels{ii,1};
    o_m='';s_m='';sh_o='';sh_s='';sh_b1='';sh_b2='';
    sh_d='';sh_tvp='';sh_vo='';sh_ssm='';
    sh_pl='';o_d='';
    is_def=equation_type(ii)==2;
    is_tvp=equation_type(ii)==3;
    is_sseq=equation_type(ii)==4;
    is_planner=equation_type(ii)==5;
    ss_is_on=false;
    
    for jj=1:size(eq_i,2)
        item=eq_i{1,jj};
        lead_or_lag=eq_i{2,jj};
        [status,pos]=parser.determine_status(dictionary,item);
        switch status
            case 'y'
                if is_sseq || is_planner
                    index=pos;
                else
                    index=abs(lead_or_lag-2);
                    index=dictionary.Lead_lag_incidence(pos,index);
                end
                if is_def
                    error([mfilename,':: definitions cannot contain variables']);
                elseif is_tvp
                    sh_tvp=[sh_tvp,'y(',sprintf('%0.0f',index),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'y(',sprintf('%0.0f',index),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'y(',sprintf('%0.0f',index),')'];
                else
                    o_m=[o_m,item];
                    if lead_or_lag
                        o_m=[o_m,'{',sprintf('%0.0f',lead_or_lag),'}'];
                    end
                    sh_o=[sh_o,'y(',sprintf('%0.0f',index),')'];
                    if ss_is_on
                        sh_vo=[sh_vo,'y(',sprintf('%0.0f',index),')'];
                    else
                        sh_vo=[sh_vo,'y(',sprintf('%0.0f',index),',:)'];
                    end
                    s_m=[s_m,item];
                    sh_s=[sh_s,'y(',sprintf('%0.0f',pos),')'];
                    this_bgp_lead=overall_max_lead_lag+lead_or_lag;
                    switch lead_or_lag
                        case -1
                            sh_b1=[sh_b1,'(','y(',sprintf('%0.0f',pos),')-y(',sprintf('%0.0f',pos+orig_endo_nbr),'))'];
                        case 0
                            sh_b1=[sh_b1,'y(',sprintf('%0.0f',pos),')'];
                        case 1
                            sh_b1=[sh_b1,'(','y(',sprintf('%0.0f',pos),')+y(',sprintf('%0.0f',pos+orig_endo_nbr),'))'];
                    end
                    switch this_bgp_lead
                        case 0
                            sh_b2=[sh_b2,'y(',sprintf('%0.0f',pos),')'];
                        case 1
                            sh_b2=[sh_b2,'(','y(',sprintf('%0.0f',pos),')+y(',sprintf('%0.0f',pos+orig_endo_nbr),'))'];
                        otherwise
                            sh_b2=[sh_b2,'(','y(',sprintf('%0.0f',pos),')+',sprintf('%0.0f',this_bgp_lead),'*y(',sprintf('%0.0f',pos+orig_endo_nbr),'))'];
                    end
                end
                ss_is_on=false;
            case 'x'
                if is_def
                    error([mfilename,':: definitions cannot contain variables']);
                elseif is_tvp
                    error([mfilename,':: endogenous switching probabilities cannot contain shocks. Use auxiliary variables if necessary']);
                elseif is_sseq
                    sh_ssm=[sh_ssm,'x(',sprintf('%0.0f',pos),')'];
                elseif is_planner
                    error([mfilename,':: planner objectives cannot include shocks. Use auxiliary variables if necessary'])
                else
                    o_m=[o_m,item];
                    sh_o=[sh_o,'x(',sprintf('%0.0f',pos),')'];
                    sh_vo=[sh_vo,'x(',sprintf('%0.0f',pos),',:)'];
                    s_m=[s_m,item];
                    sh_s=[sh_s,'x(',sprintf('%0.0f',pos),')'];
                    sh_b1=[sh_b1,'x(',sprintf('%0.0f',pos),')'];
                    sh_b2=[sh_b2,'x(',sprintf('%0.0f',pos),')'];
                end
            case 'param'
                if is_def
                    sh_d=[sh_d,'param(',sprintf('%0.0f',pos),')'];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,'param(',sprintf('%0.0f',pos),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'param(',sprintf('%0.0f',pos),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'param(',sprintf('%0.0f',pos),')'];
                else
                    o_m=[o_m,item];
                    sh_o=[sh_o,'param(',sprintf('%0.0f',pos),')'];
                    sh_vo=[sh_vo,'param(',sprintf('%0.0f',pos),')'];
                    s_m=[s_m,item];
                    sh_s=[sh_s,'param(',sprintf('%0.0f',pos),')'];
                    sh_b1=[sh_b1,'param(',sprintf('%0.0f',pos),')'];
                    sh_b2=[sh_b2,'param(',sprintf('%0.0f',pos),')'];
                end
            case 'def'
                if is_def
                    sh_d=[sh_d,'def(',sprintf('%0.0f',pos),')'];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,'def(',sprintf('%0.0f',pos),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'def(',sprintf('%0.0f',pos),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'def(',sprintf('%0.0f',pos),')'];
                else
                    o_m=[o_m,item];
                    s_m=[s_m,item];
                    sh_o=[sh_o,'def(',sprintf('%0.0f',pos),')'];
                    sh_vo=[sh_vo,'def(',sprintf('%0.0f',pos),')'];
                    sh_s=[sh_s,'def(',sprintf('%0.0f',pos),')'];
                    sh_b1=[sh_b1,'def(',sprintf('%0.0f',pos),')'];
                    sh_b2=[sh_b2,'def(',sprintf('%0.0f',pos),')'];
                end
            case 'tvp'
                if ~is_tvp
                    error([mfilename,':: endogenous probability cannot be used in the model equations'])
                end
                sh_tvp=[sh_tvp,item];
            otherwise
                if strcmp(item,'steady_state')
                    % then two blocks later is the name of the variable
                    Vss=eq_i{1,jj+2};
                    pos=find(strcmp(Vss,{dictionary.orig_endogenous.name}));
                    eq_i{1,jj+2}=sprintf('%0.0f',pos);
                end
                if is_def
                    sh_d=[sh_d,item];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,item];
                elseif is_sseq
                    sh_ssm=[sh_ssm,item];
                elseif is_planner
                    sh_pl=[sh_pl,item];
                else
                    o_m=[o_m,item];
                    s_m=[s_m,item];
                    sh_o=[sh_o,item];
                    sh_vo=[sh_vo,item];
                    sh_s=[sh_s,item];
                    sh_b1=[sh_b1,item];
                    sh_b2=[sh_b2,item];
                end
        end
    end
    if is_def
        shadow_definitions=[shadow_definitions;{sh_d}]; % put back the semicolon as this is going to be evaluated
        orig_definitions=[orig_definitions;{o_d}];
    elseif is_tvp
        shadow_tvp=[shadow_tvp;{sh_tvp}];
    elseif is_sseq
        static.steady_state_shadow_model=[static.steady_state_shadow_model;{sh_ssm}];
        % put back the semicolon as this is going to be evaluated
    elseif is_planner
        dictionary.planner_system.shadow_model=[dictionary.planner_system.shadow_model;{sh_pl}];
    else
        dynamic.model=[dynamic.model;{o_m}];
        static.model=[static.model;{s_m}];
        dynamic.shadow_model=[dynamic.shadow_model;{parser.replace_steady_state_call(sh_o)}];
        static.shadow_model=[static.shadow_model;{parser.replace_steady_state_call(sh_s)}];
        static.shadow_BGP_model=[static.shadow_BGP_model;{parser.replace_steady_state_call(sh_b1)};{parser.replace_steady_state_call(sh_b2)}];
    end
end

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
%% load transition matrices and transorm markov chains
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
    dictionary.shadow_transition_matrix,...
    dictionary.markov_chains,...
    myifelseif...
    ]=parser.transition_probabilities(...
    dictionary.input_list,...
    dictionary.parameters,dictionary.markov_chains,...
    shadow_tvp,bug_fix_shadow_defs,probability_of_commitment);
% flag for endogenous probabilities
%---------------------------------
dictionary.is_endogenous_switching_model=any(dictionary.markov_chains.chain_is_endogenous);
%% symbolic forms
%
disp(' ')
disp('Now computing symbolic derivatives...')
disp(' ')
orig_exo_nbr=numel(dictionary.exogenous);
param_nbr = numel(dictionary.parameters);
switching_ones=find([dictionary.parameters.is_switching]);
% dynamic model wrt y and x
%--------------------------
wrt={'y',1:nnz(dictionary.Lead_lag_incidence)
    'x',1:orig_exo_nbr
    'param',switching_ones}; % <-- switching parameters only
myIncidence=logical([dictionary.Lead_lag_incidence(:)',1:orig_exo_nbr,1:numel(switching_ones)]);
myPartitions={'p','c','m','e','t';
    orig_endo_nbr,orig_endo_nbr,orig_endo_nbr,orig_exo_nbr,numel(switching_ones)
    };
%     [model_derivatives,~,zeroth_order,numEqtns,numVars,jac_toc]=differentiate_system(...
% % profile off, profile on
[model_derivatives,numEqtns,numVars,jac_toc]=differentiate_system(...
    parser.burry_probabilities(dynamic.shadow_model,myifelseif),... % 
    dictionary.input_list,myIncidence,wrt,myPartitions,2);
disp([mfilename,':: 1st and 2nd-order derivatives of dynamic model wrt y(+0-), x and theta ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])
% % profile off, profile viewer
% % keyboard

% static model wrt y
%--------------------
wrt={'y',1:orig_endo_nbr};
myIncidence=logical(1:orig_endo_nbr)';
myPartitions={'c';orig_endo_nbr};
[static_model_derivatives,numEqtns,numVars,jac_toc]=differentiate_system(...
    static.shadow_model,...
    dictionary.input_list,myIncidence,wrt,myPartitions,1);
disp([mfilename,':: 1st-order derivatives of static model wrt y(0). ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])

% balanced growth path model wrt y
%---------------------------------
wrt={'y',1:2*orig_endo_nbr};
myIncidence=logical(1:2*orig_endo_nbr)';
myPartitions={'c';2*orig_endo_nbr};
[static_bgp_model_derivatives,numEqtns,numVars,jac_toc]=differentiate_system(...
    static.shadow_BGP_model,...
    dictionary.input_list,myIncidence,wrt,myPartitions,1);
disp([mfilename,':: 1st-order derivatives of static BGP model wrt y(0). ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])

% dynamic model wrt param
%------------------------
wrt={'param',1:param_nbr};
myIncidence=logical(1:param_nbr)';
myPartitions={'c';param_nbr};
ppdd=@(x)x;%dynamic.shadow_model;
if DefaultOptions.definitions_in_param_differentiation
    ppdd=@(x)parser.replace_definitions(x,shadow_definitions);
end
[param_derivatives,numEqtns,numVars,jac_toc]=differentiate_system(...
    ppdd(dynamic.shadow_model),...
    dictionary.input_list,myIncidence,wrt,myPartitions,1);
disp([mfilename,':: first-order derivatives of dynamic model wrt param. ',...
    sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.4f',jac_toc),' seconds'])
%% push derivatives in to functions
% for the moment, put all in the same matrix
dynamic.model_derivatives=struct(...
    'Endogenous_Shocks_Parameters',{model_derivatives},...
    'Parameters',{param_derivatives},...
    'StaticEndogenous',{static_model_derivatives},...
    'Static_BGP_Endogenous',{static_bgp_model_derivatives}...
    );

if is_model_with_planner_objective
    wrt={'y',1:size(dictionary.Lead_lag_incidence,1)};
    myIncidence=logical(1:orig_endo_nbr)';
    myPartitions={'c';orig_endo_nbr};
    % do not chop the output because we are going to augment the function
    % with the loss, the degree of commitment and the discount. We achieve
    % this by setting the last argument of the function to false.
    [LossComDiscHessJac,numEqtns,numVars,jac_toc]=differentiate_system(...
        dictionary.planner_system.shadow_model(1),...
        dictionary.input_list,myIncidence,wrt,myPartitions,2);
    disp([mfilename,':: 1st and 2nd-order derivatives of planner objective wrt y(0). ',...
        sprintf('%0.0f',numEqtns),' equations and ',sprintf('%0.0f',numVars),' variables :',sprintf('%0.0f',jac_toc),' seconds'])
    % add the loss, the commitment degree and discount
    %-------------------------------------------------
    shadow_model=dictionary.planner_system.shadow_model;
    tmp=code2func(strrep(strrep(shadow_model,'discount-',''),'commitment-',''),dictionary.input_list,true);
    % add some empty fields to conform with the differentiation
    tmp.map=[]; tmp.partitions=[];tmp.maxcols=[];
    % put all in a single structure
    LossComDiscHessJac=[tmp,LossComDiscHessJac];
    dictionary.planner_system.LossComDiscHessJac=LossComDiscHessJac;
end
%% Add final variables list to the dictionary as they may differ earlier ones
unsorted_endogenous=dictionary.orig_endogenous;

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

if dictionary.is_optimal_policy_model
    if dictionary.is_sticky_information_model
        error([mfilename,':: you are not allowed to solve a sticky information model with loose commitment'])
    end
    if dictionary.is_hybrid_expectations_model
        error([mfilename,':: you are not allowed to solve a hybrid expectations model with loose commitment'])
    end
    
    for eq=1:dictionary.NumberOfEquations
        new_var=struct('name',['mult_',sprintf('%0.0f',eq)],'tex_name','',...
            'max_lead',-equations_maxLag_maxLead(eq,1),... % lags govern the leads
            'max_lag',-equations_maxLag_maxLead(eq,2),... % keads govern the lags
            'is_log_var',false,'is_auxiliary',false);
        for ifield=1:numel(useless_fields)
            new_var.(useless_fields{ifield})=nan;
        end
        unsorted_endogenous=[unsorted_endogenous,new_var];
    end
    added=[-equations_maxLag_maxLead(:,1),...
        ones(dictionary.NumberOfEquations,1),...
        -equations_maxLag_maxLead(:,2)];
logical_incidence=[dictionary.Lead_lag_incidence;added];
else
    added=0;
    assert(numel(dictionary.orig_endogenous)==sum(equation_type==1),...
        '# equations different from # endogenous variables')
    if dictionary.is_sticky_information_model
        dictionary.forward_looking_ids=find(dictionary.Lead_lag_incidence(:,3));
        for ii=1:numel(dictionary.forward_looking_ids)
            id=dictionary.forward_looking_ids(ii);
            new_var=struct('name',['SI_',dictionary.orig_endogenous{id}],'tex_name','',...
                'max_lead',0,'max_lag',0,...
                'is_log_var',false,'is_auxiliary',true);
            for ifield=1:numel(useless_fields)
                new_var.(useless_fields{ifield})=nan;
            end
            unsorted_endogenous=[unsorted_endogenous,new_var];
        end
        added=numel(dictionary.forward_looking_ids);
    end
    % this will be used for the determination of the status of the variables.
    % the added variables are given status of static, which is misleading when
    % added>0
    logical_incidence=[dictionary.Lead_lag_incidence;repmat([0,1,0],added,1)];
end
% now we can resort the final variables
[~,tags]=sort({unsorted_endogenous.name});
logical_incidence=logical_incidence(tags,:);
dictionary.endogenous=unsorted_endogenous(tags);
% this will be used to reorder the solution of loose commitment or sticky
% information
dictionary.reordering_index=locate_variables({dictionary.endogenous.name},{unsorted_endogenous.name});
dictionary.reordering_index=transpose(dictionary.reordering_index(:));
clear unsorted_endogenous

%% give greek names to endogenous, exogenous, parameters
dictionary.endogenous=parser.greekify(dictionary.endogenous);
dictionary.orig_endogenous=parser.greekify(dictionary.orig_endogenous);
dictionary.exogenous=parser.greekify(dictionary.exogenous);
dictionary.parameters=parser.greekify(dictionary.parameters);

dictionary.is_purely_forward_looking_model=false;
dictionary.is_purely_backward_looking_model=false;
dictionary.is_hybrid_model=any(dictionary.Lead_lag_incidence(:,1)) && any(dictionary.Lead_lag_incidence(:,3));
if ~dictionary.is_hybrid_model
    if any(dictionary.Lead_lag_incidence(:,1))
        dictionary.is_purely_forward_looking_model=true;
    elseif any(dictionary.Lead_lag_incidence(:,3))
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
dictionary.endogenous.is_frwrd_looking=sparse((logical_incidence(:,1) & ~logical_incidence(:,3))');
dictionary.endogenous.is_log_var=sparse([endogenous.is_log_var]);
dictionary.endogenous.is_auxiliary=sparse([endogenous.is_auxiliary]);
clear endogenous logical_incidence

exogenous=dictionary.exogenous;
dictionary.exogenous=struct();
dictionary.exogenous.name={exogenous.name};
dictionary.exogenous.tex_name={exogenous.tex_name};
dictionary.exogenous.is_observed=sparse(ismember(dictionary.exogenous.name,{dictionary.observables.name}));
dictionary.exogenous.number=full([sum(~dictionary.exogenous.is_observed),sum(dictionary.exogenous.is_observed)]);
dictionary.exogenous.is_in_use=sparse([exogenous.is_in_use]);
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

dictionary.model_derivatives=dynamic.model_derivatives;
dictionary.steady_state_shadow_model=static.steady_state_shadow_model;
%  char({dictionary.equations.shadow_balanced_growth_path})
definitions=dictionary.definitions;
dictionary.definitions=struct();
dictionary.definitions.dynamic={definitions.model};
dictionary.definitions.dynamic=dictionary.definitions.dynamic(:);
dictionary.definitions.shadow_dynamic={definitions.shadow};
dictionary.definitions.shadow_dynamic=dictionary.definitions.shadow_dynamic(:);
dictionary.definitions.number=numel(definitions);
clear definitions static dynamic

dictionary=orderfields(dictionary);

end

function [derivs,numEqtns,numVars,jac_toc]=differentiate_system(myfunc,input_list,incidence,wrt,Partitions,order)

with_respect_to=extend_differentiation_list(incidence,wrt);

myfunc=analytical_symbolic_form(myfunc,input_list,'symbolic');

% list of symbols
symb_list=collect_symbolic_list(myfunc,strcat(input_list,'_'));
% force 's0' and 's1' to enter the list
state_inputs={'s0','s1'};
input_list=input_list(:)';
for ii=1:numel(state_inputs)
    if ~any(strcmp(symb_list,state_inputs{ii}))
        symb_list=[symb_list,state_inputs{ii}];
    end
    if ~any(strcmp(input_list,state_inputs{ii}))
        input_list=[input_list,state_inputs{ii}];
    end
end
% sorting will be useful if we need to debug
symb_list=sort(symb_list);

args=planar.initialize(symb_list,with_respect_to);
numEqtns=numel(myfunc);
numVars=numel(with_respect_to);
for ifunc=1:numEqtns
    [occur,myfunc{ifunc}]=find_occurrences(myfunc{ifunc},symb_list);
    % re-create the function
    var_occur=symb_list(occur);
    argfun=cell2mat(strcat(var_occur,','));
    myfunc{ifunc}=str2func(['@(',argfun(1:end-1),')',myfunc{ifunc}]);
    arg_occur=args(occur);
    myfunc{ifunc}=myfunc{ifunc}(arg_occur{:});
end
verbose=true;
compact_derivatives=true;
tic
derivs=planar.differentiate(myfunc,order,Partitions,verbose);
for oo=1:order
    if verbose
        tic
    end
    derivs(oo)=planar.print(derivs(oo),false);
    if verbose
        fprintf(1,'printing of derivatives at order %0.0f done in %0.4f seconds\n',oo,toc);
        tic
    end
    derivs(oo)=planar.derivatives2functions(derivs(oo),input_list,compact_derivatives);
    if verbose
        fprintf(1,'derivatives to functions at order %0.0f done in %0.4f seconds\n',oo,toc);
    end
end
jac_toc=toc;

myderivs=struct();
ff=fieldnames(derivs);
for oo=1:order
    for ifield=1:numel(ff)
        f1=ff{ifield};
        if strcmp(f1,'derivatives')
            f1='functions';
        end
        myderivs(oo).(f1)=derivs(oo).(ff{ifield});
    end
end
derivs=myderivs;

    function wrt=extend_differentiation_list(incidence,wrt0)
        
        set_differentiation_list();
        
        incidence=incidence(:);
        n=numel(incidence);
        wrt=cell(1,n);
        wrt(incidence)=wrt0;
        iter=0;
        for ii_=1:n
            if isempty(wrt{ii_})
                iter=iter+1;
                wrt{ii_}=['zZzZzZz_',sprintf('%0.0f',iter)];
            end
        end
        
        function set_differentiation_list()
            reprocess=size(wrt0,2)==2 && isnumeric(wrt0{1,2});
            if reprocess
                with_respect_to=cell(1,300);
                jter=0;
                for iii=1:size(wrt0,1)
                    digits=wrt0{iii,2};
                    xx=wrt0{iii,1};
                    for id=1:numel(digits)
                        jter=jter+1;
                        if jter==size(with_respect_to,2)
                            with_respect_to{end+300}={};
                        end
                        with_respect_to{jter}=[xx,'_',sprintf('%0.0f',digits(id))];
                    end
                end
                wrt0=with_respect_to(1:jter);
            end
        end
    end
end
