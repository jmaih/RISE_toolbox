function dictionary=CompileModelFile(FileName,varargin)
% by the way, I can still declare exogenous and make them observable at the
% same time. The exogenous that are observed are determisitic. This opens
% the door for estimating partial equilibrium models

DefaultOptions=struct('definitions_in_param_differentiation',true);
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
equation=cell(2,0);
tex_name='';
dollar_active=false;
last_block_id=[];
function_on=0;
time_on=false;
last_status='';
fill_time='';
time_opening='';
def_flag=false;
endo_switch_flag=false;

dictionary=struct();
dictionary.is_linear_model=false;
static.is_imposed_steady_state=false;
static.is_unique_steady_state=false;
% those names will be used in functions and in this specific order
dictionary.input_list={'y','x','param','ss','def'};

%% set various blocks

DELIMITERS=[char([9:13,32]),'[]{}(),;=+-*/^@><'];
% first output: the dictionary.filename
is_svar_model= isstruct(FileName);
if is_svar_model
    svar_fields={'model','var','varexo'}; 
    RawFile=[];
    % keeping the same keywords as the ones that would be used in a model
    % file
    iter_svar=0;
    for ifield=1:numel(svar_fields)
        thisField=svar_fields{1,ifield};
        if ~isfield(FileName,thisField)
            if strcmp(thisField,'varexo')
                continue
            else
                error([thisField,' must be a field when declaring a svar'])
            end
        end
        if strcmp(thisField,'model') && ~isequal(FileName.(thisField),'svar')
            error('The entry for model should be ''svar'' when declaring a svar')
        end
        if ismember(thisField,{'var','varexo'})
            iter_svar=iter_svar+1;
            if ischar(FileName.(thisField))
                FileName.(thisField)=cellstr(FileName.(thisField));
            end
            tmp=FileName.(thisField);
            tmp=tmp(:)';
            tmp=strcat(tmp,',');
            tmp=cell2mat(tmp);
            % create a line for the endogenous
            RawFile=[RawFile
                {[thisField,',',tmp(1:end-1),';'],'svar',iter_svar}];
        end
    end
    dictionary.filename='svar';
else
    FileName(isspace(FileName))=[];
    loc=strfind(FileName,'.');
    if isempty(loc)
        dictionary.filename=FileName;
        if exist([FileName,'.dyn'],'file')
            FileName=[FileName,'.dyn'];
        elseif exist([FileName,'.mod'],'file')
            FileName=[FileName,'.mod'];
        else
            error([mfilename,':: ',FileName,'.mod or ',FileName,'.dyn not found'])
        end
    else
        ext=FileName(loc:end);
        if strcmp(ext,'.dyn')||strcmp(ext,'.mod')
        else
            error([mfilename,':: Input file is expected to be a .dyn or .mod file'])
        end
        dictionary.filename=FileName(1:loc-1);
    end

% read file and remove comments
RawFile=read_file(FileName);
end

blocks=struct('name',{'log_vars','orig_endogenous','exogenous','parameters',...
    'observables','model','steady_state_model','parameterization',...
    'planner_objective','exogenous_definition','parameter_restrictions'},...
    'trigger',{'log_vars','var','varexo','parameters','varobs','model',...
    'steady_state_model','parameterization','planner_objective',...
    'exogenous_definition','parameter_restrictions'},...
    'active',num2cell(false(1,11)),...
    'listing',{cell(0,2),cell(0,2),cell(0,2),cell(0,2),cell(0,2),cell(0,3),...
    cell(0,3),cell(0,3),cell(0,3),cell(0,3),cell(0,3)});

% 'orig_endogenous','exogenous','parameters','observables' blocks are just
% declarations. The corresponding blocks are to hold name and tex_name
% respectively

% construct cells with 3 columns that will hold line_number, equation,
% dictionary.filename respectively


NumberOfLines=size(RawFile,1);
iline=0;
while iline<NumberOfLines
    iline=iline+1;
    rawline=RawFile{iline,1};
    file_name=RawFile{iline,2};
    line_number=RawFile{iline,3};
    
    [tok,rest]=strtok(rawline,DELIMITERS);
    [is_trigger,blocks,last_block_id]=check_block(tok,blocks,last_block_id,file_name,line_number);
    if isempty(last_block_id)
        error([mfilename,':: string ''',tok,''' in ''',file_name,''' at line ',int2str(line_number),' does not belong to a known block '])
    end
    current_block_name=blocks(last_block_id).name;
    switch current_block_name
        case {'log_vars','orig_endogenous','exogenous','parameters','observables'}
            blocks(last_block_id)=construct_list(blocks(last_block_id),rawline,tok);
        case {'model','steady_state_model','parameterization',...
                'planner_objective','exogenous_definition','parameter_restrictions'}
            if is_trigger && ismember(current_block_name,...
                    {'parameterization','planner_objective','parameter_restrictions'})
                % remove the trigger and parse the rest
                rawline=rest;
                rawline(isspace(rawline))=[];
                if strcmp(rawline,';')
                    rawline=[];
                end
            end
            if isempty(rawline)
                continue
            end
            end_game=(logical(strcmp(tok,'end')) && ismember(current_block_name,...
                {'model','steady_state_model','parameterization',...
                'exogenous_definition','parameter_restrictions'}))||...
                (~isempty(strfind(rawline,';')) && strcmp(current_block_name,'planner_objective'));
            if end_game
                blocks(last_block_id).active=false;
                if strcmp(current_block_name,'planner_objective')
                    blocks(last_block_id).listing=[blocks(last_block_id).listing;{line_number,rawline,file_name}];
                end
            else
                blocks(last_block_id).listing=[blocks(last_block_id).listing;{line_number,rawline,file_name}];
            end
    end
end
if any([blocks.active])
    disp({blocks([blocks.active]).name})
    error([mfilename,':: block ',blocks([blocks.active]).name,' is not closed'])
end

if is_svar_model
    % declare the shocks and the observables as well
    blknames={blocks.name};
    exoblock=strcmp(blknames,'exogenous');
    varobsblock=strcmp(blknames,'observables');
    endoblock=strcmp(blknames,'orig_endogenous');
    exoListing=blocks(endoblock).listing;
    varobsListing=blocks(endoblock).listing;
    for ilist=1:size(exoListing,1)
        exoListing{ilist,1}=['EPS_',exoListing{ilist,1}];
        if ~isempty(exoListing{ilist,2})
            exoListing{ilist,2}=[exoListing{ilist,2},' shock'];
            varobsListing{ilist,2}='';
        end
    end
    observed_exo=blocks(exoblock).listing;
    blocks(exoblock).listing=[observed_exo;exoListing];
    blocks(varobsblock).listing=[observed_exo;varobsListing];
end
clear current_block_name
%% Populate the dictionary
dictionary.definitions={};
% unlike the declared list, definitions will never be sorted
dictionary.known_words={'steady_state','argzero','x0_','x1_','param_obj','commitment','discount',...
    'log','exp','cos','sin','normpdf','normcdf'};
dictionary.known_words=[dictionary.known_words,strcat({'xx_ssmdef_'},cellstr(int2str((1:9)'))')];
% dictionary.steady_state_parameters={};
dictionary.time_varying_probabilities={};
dictionary.symbols={'#','!',')','(','}','{',']','[',',',';','.','=','@'};
dictionary.add_operators={'+','-'};
dictionary.mult_operators={'*','^','/'};
dictionary.relational_operators={'<','>'};
distr_list=what(['classes',filesep,'+distributions']);
for ilist=1:numel(distr_list)
    if ~isempty(distr_list(ilist).m)
        mydistrlist=distr_list(ilist).m;
        break
    end
end
% for some reason I do not understand, this sometimes returns a 2 x 1
% structure instead of a 1 x 1. I have experienced it when using parallel
% computing, but not otherwise.
dictionary.Distributions=strrep(mydistrlist,'.m','_pdf');
% % % dictionary.Distributions={'uniform_pdf','unif_pdf','normal_pdf','norm_pdf',...
% % %     'gamma_pdf','gam_pdf','beta_pdf','invg_pdf','inv_gamma_pdf',...
% % %     'inv_gamma2_pdf','dirichlet_pdf'};
dictionary.syntax_special={'y]','param]',')]','x]','[(',']=','@f','(@','])',...
    '([','[n','n]',...
    '>=','<=','param>','>param',')>','>(','>f','(cn','cn)','cn,',',cn',',n',...
    '>n','<n'};
% last line added for the parameter restrictions block
dictionary.syntax_function={'y,','x,','param,','),','f,','n,',',y',',x',...
    ',param',',f',',n',',('};
dictionary.syntax_time={'y(','y{','x(','x{','n}','{+','}+','}*','{n','})','};','+}','}='};
dictionary.syntax_typical={'y)','y+','y*','y;','x)','x+','x*','x;','param)','def)',...
    'param+','def+','param*','def*','param;','def;','))',')+',')*',');','(y','(x',...
    '(param','(def','((','(+','(f','(n','+y','+x','+param','+def','+(','+f','+n',...
    '*y','*x','*param','*def','*f','*n','*(','f)','f(','n.','.n','f;','n)','n+',...
    'n*','nn','n;','=.','=param','=def','=x','=y','=f','=(','=+','=n',...
    'param=','def=','x=','y=','f=','n=',')=','tvp='};
% endogenous(y),exogenous(x),parameter(param), known word or function (f),
% number(n), definition(def)

% add the list of endogenous,exogenous,parameters and observables and
% remove them from the blocks
dic_items={'orig_endogenous','exogenous','parameters','observables','log_vars'};
for ii=1:numel(dic_items)
    loc=strcmp(dic_items{ii},{blocks.name});
    dictionary.(dic_items{ii})=struct('name',{},'tex_name',{});
    % sort variables and parameters right here right now
    [~,tags]=sort(blocks(loc).listing(:,1));
    for jj=1:numel(tags)
        dictionary.(dic_items{ii})(jj).name=blocks(loc).listing{tags(jj),1};
        dictionary.(dic_items{ii})(jj).tex_name=blocks(loc).listing{tags(jj),2};
        if ~ismember(dic_items{ii},{'observables','log_vars'})
            % parameters can have leads in the framework of Rubio-Ramirez,
            % Waggoner and Zha (2012). While my framework takes it as given
            % and infers the leads as the parameters attached to the lead
            % terms, I make provision for the RRWZ case in case I happen to
            % change my mind on that some day.
            dictionary.(dic_items{ii})(jj).max_lead=0;
            dictionary.(dic_items{ii})(jj).max_lag=0;
        end
        if strcmp(dic_items{ii},'parameters')
            dictionary.(dic_items{ii})(jj).is_switching=false;
            dictionary.(dic_items{ii})(jj).is_measurement_error=false;
        end
    end
    % remove item from block
    blocks(loc)=[];
end
for ii=1:numel(dic_items)
    vi={dictionary.(dic_items{ii}).name};
    if strcmp(dic_items{ii},'observables')
        % check that all observables are either endogenous or exogenous
        badguys=~ismember(vi,[{dictionary.orig_endogenous.name},{dictionary.exogenous.name}]);
        if any(badguys)
            disp(vi(badguys))
            error([mfilename,':: the observable variables above are not declared either as exogenous or endogenous'])
        end
    end
    if strcmp(dic_items{ii},'log_vars')
        % check that all log_vars are endogenous
        badguys=~ismember(vi,{dictionary.orig_endogenous.name});
        if any(badguys)
            disp(vj(badguys))
            error([mfilename,':: the log variables above are not declared as endogenous'])
        end
    end
    % check that the same name is not found in other places
    for jj=ii+1:numel(dic_items)
        if (~ismember('observables',{dic_items{ii},dic_items{jj}}) && ~ismember('log_vars',{dic_items{ii},dic_items{jj}}))||...
                (ismember('observables',{dic_items{ii},dic_items{jj}}) && ismember('log_vars',{dic_items{ii},dic_items{jj}}))
            vj={dictionary.(dic_items{jj}).name};
            badguys=ismember(vj,vi);
            if any(badguys)
                disp(vj(badguys))
                more_string='';
                if (ismember('observables',{dic_items{ii},dic_items{jj}}) && ismember('log_vars',{dic_items{ii},dic_items{jj}}))
                    more_string=' Use auxiliary variables for measurement if necessary';
                end
                error([mfilename,':: atoms above cannot be both ''',dic_items{ii},''' and ''',dic_items{jj},'''. ',more_string])
            end
        end
    end
end
%% Markov chain information (constant parameters)
% Initialize markov chain information with constant parameters
dictionary.MarkovChains=update_markov_chains_info;

% update markov chain info by ripping through the parameter list. this does
% not depend on anything else than the parameter list. Later on, we will
% update this information with time-varying probabilities.
for ii=1:numel(dictionary.parameters)
    dictionary.MarkovChains=update_markov_chains_info(dictionary.MarkovChains,dictionary.parameters(ii).name);
end

%% add the chain names to the dictionary, for parsing purposes
dictionary.chain_names=dictionary.MarkovChains(1,:);

%% Model block
% now with the endogenous, exogenous, parameters in hand, we can process
% the model block

% it may be useful to keep track of the parameters and shocks that are
% actually in use
is_in_use_shock=false(numel(dictionary.exogenous),1);
is_in_use_parameter=false(numel(dictionary.parameters),1);

current_block_id=find(strcmp('model',{blocks.name}));
Model_block=cell(0,1);
for ii=1:size(blocks(current_block_id).listing,1)
    Model_block=capture_equations(Model_block,blocks(current_block_id).listing(ii,:),'model');
end
if ~ is_svar_model && isempty(Model_block)
    error([mfilename,':: no model declared'])
end
% remove item from block
blocks(current_block_id)=[];

%% add more equations if necessary, such that the max_lag and max_lead are 1
% replace exogenous lags with auxiliary endogenous variables
% update the number of endogenous variables accordingly
% Also keep the same equations to be added to the steady state model
auxiliary_steady_state_equations=cell(0,1);
orig_endogenous_current=dictionary.orig_endogenous; % these are the variables without the augmentation add-on
for ii=1:numel(dictionary.orig_endogenous)
    if dictionary.orig_endogenous(ii).max_lead<2
        continue
    end
    vname=dictionary.orig_endogenous(ii).name;
    lead_i=dictionary.orig_endogenous(ii).max_lead;
    vold=vname;
    for i2=2:lead_i
        new_var=struct('name',[vname,'_AUX_F_',int2str(i2-1)],'tex_name','','max_lead',1,'max_lag',0);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vold,1}',{';',[]}']}]; %
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{vold,0}',{';',[]}']}]; %
        if i2-1==1 % set the lead to 1
            dictionary.orig_endogenous(ii).max_lead=1;
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
        new_var=struct('name',[vname,'_0'],'tex_name','','max_lead',0,'max_lag',-1);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vname,0}',{';',[]}']}]; %
        % The steady state is computed with zero shocks. and so instead of
        % setting the name of the shock, I write 0
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{'0',0}',{';',[]}']}]; %
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
        new_var=struct('name',[vname,'_AUX_L_',int2str(i2-1)],'tex_name','','max_lead',0,'max_lag',-1);
        Model_block=[Model_block;
            {[{new_var.name,0}',{'-',[]}',{vold,-1}',{';',[]}']}]; %
        auxiliary_steady_state_equations=[auxiliary_steady_state_equations
            {[{new_var.name,0}',{'=',[]}',{vold,0}',{';',[]}']}]; %
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

% sorting the endogenous switching probabilities is more or less useless
dictionary.time_varying_probabilities=sort(dictionary.time_varying_probabilities);

%% Markov chain information (time-varying probabilities)
% update the endogenous switching variables, but now setting the
% flag to true
for ii=1:numel(dictionary.time_varying_probabilities)
    dictionary.MarkovChains=update_markov_chains_info(dictionary.MarkovChains,dictionary.time_varying_probabilities{ii},true);
end

%% update the chain names in the dictionary, for parsing purposes
% this will help parse, for instance, the parameter restrictions block
dictionary.chain_names=dictionary.MarkovChains(1,:);

%% Now we re-write the model and update leads and lags
% at the same time, construct the incidence and occurrence matrices
orig_endo_nbr=numel(dictionary.orig_endogenous);

number_of_equations=numel(Model_block);
Occurrence=false(number_of_equations,orig_endo_nbr,3);
equation_type=ones(number_of_equations,1);
for ii=1:number_of_equations
    eq_i= Model_block{ii};
    if ismember(eq_i{1,1},dictionary.definitions) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=2;
    elseif ismember(eq_i{1,1},dictionary.time_varying_probabilities) && strcmp(eq_i{1,2}(1),'=')
        equation_type(ii)=3;
    end
    for i2=1:size(eq_i,2)
        if ~isempty(eq_i{2,i2})
            vname=eq_i{1,i2};
            status=determine_status(vname);
            time=-1; new_item=false;
            if strcmp(status,'x')&& abs(eq_i{2,i2})>0
                new_item=true;
                if abs(eq_i{2,i2})==1
                    % change the name, but not the lag structure
                    eq_i{1,i2}=[vname,'_0'];
                else
                    % change both the name and the lag structure
                    eq_i{1,i2}=[vname,'_0','_AUX_L_',int2str(abs(eq_i{2,i2})-1)];
                end
            elseif strcmp(status,'y')&& abs(eq_i{2,i2})>1
                new_item=true;
                % change both the name and the lag structure
                if eq_i{2,i2}>0
                    eq_i{1,i2}=[vname,'_AUX_F_',int2str(abs(eq_i{2,i2})-1)];
                    time=1;
                elseif eq_i{2,i2}<0
                    eq_i{1,i2}=[vname,'_AUX_L_',int2str(abs(eq_i{2,i2})-1)];
                end
            end
            if new_item
                eq_i{2,i2}=time;
            end
            var_loc=strcmp(eq_i{1,i2},{dictionary.orig_endogenous.name});
            lag_or_lead=eq_i{2,i2}+2;
            if equation_type(ii)==2
                error([mfilename,':: equation (',int2str(ii),') detected to be a definition cannot contain variables'])
            elseif equation_type(ii)==3 && ismember(lag_or_lead,[3,1])
                error([mfilename,':: equation (',int2str(ii),') detected to describe endogenous switching cannot contain leads or lags'])
            end
            Occurrence(ii,var_loc,lag_or_lead)=true;
        end
    end
    Model_block{ii}= eq_i;
end
% keep only the structural equations
Occurrence=Occurrence(equation_type==1,:,:);

%% Steady state Model block
% now with the endogenous, exogenous, parameters in hand, we can process
% the steady state model block

current_block_id=find(strcmp('steady_state_model',{blocks.name}));
SteadyStateModel_block=cell(0,1);
for ii=1:size(blocks(current_block_id).listing,1)
    SteadyStateModel_block=capture_equations(SteadyStateModel_block,blocks(current_block_id).listing(ii,:),'steady_state_model');
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
ExogenousDefinition_block=cell(0,1);
for ii=1:size(blocks(current_block_id).listing,1)
    ExogenousDefinition_block=capture_equations(ExogenousDefinition_block,blocks(current_block_id).listing(ii,:),'exogenous_definition');
end
% remove item from block
blocks(current_block_id)=[];
% the equations have been validated, now rebuild them and keep a list of
% the variables defined
DefinedExoList=cell(1,0);%{}
for ii=1:numel(ExogenousDefinition_block)
    DefinedExoList=[DefinedExoList,ExogenousDefinition_block{ii}(1,1)];
    eq_i='';
    for jj=1:size(ExogenousDefinition_block{ii},2)
        eq_i=[eq_i,ExogenousDefinition_block{ii}{1,jj}];
        if ~isempty(ExogenousDefinition_block{ii}{2,jj}) && ExogenousDefinition_block{ii}{2,jj}~=0
            eq_i=[eq_i,'{',int2str(ExogenousDefinition_block{ii}{2,jj}),'}'];
        end
    end
    ExogenousDefinition_block{ii}=eq_i;
end
% assign information to dictionary
dictionary.exogenous_equations=struct('name',DefinedExoList,'equation',transpose(ExogenousDefinition_block));
clear ExogenousDefinition_block DefinedExoList
%% optimal policy block
current_block_id=find(strcmp('planner_objective',{blocks.name}));
is_model_with_planner_objective=~isempty(blocks(current_block_id).listing);
PlannerObjective_block=cell(0,1);

dictionary.planner=[];
if is_model_with_planner_objective
    % there can't be multiple planner_objective blocks
    Commitment='commitment=1;';
    Discount='discount=.99;';
    Objective='';
    % pull everything on one line
    fichier_name_=blocks(current_block_id).listing{1,3};
    line_number_=blocks(current_block_id).listing{1,1};
    OneLiner='';
    for ii=1:size(blocks(current_block_id).listing,1)
        OneLiner=[OneLiner,blocks(current_block_id).listing{ii,2}];
        if ~ismember(blocks(current_block_id).listing{ii,1},line_number_)
            line_number_=[line_number_,blocks(current_block_id).listing{ii,1}];
        end
    end
    OneLiner(isspace(OneLiner))=[];
    discount_flag=false;
    commitment_flag=false;
    while ~isempty(OneLiner)
        [tok,rest]=strtok(OneLiner,'={},;');
        if strcmp(tok,'commitment')
            commitment_flag=true;
        elseif strcmp(tok,'discount')
            discount_flag=true;
        end
        if commitment_flag||discount_flag
            if ~strcmp(rest(1),'=')
                error([mfilename,':: expected an equality after either commitment or discount in ''',fichier_name_,''' at line ',line_number_])
            end
            rest=rest(2:end);
            % take anything before the next , or } or ;
            [tok,rest]=strtok(rest,'},;');
            if commitment_flag
                Commitment=['commitment=',tok,';'];
                commitment_flag=false;
            elseif discount_flag
                Discount=['discount=',tok,';'];
                discount_flag=false;
            end
            if strcmp(rest(1),'}')
                rest=rest(2:end);
            end
        else
            [look_behind]=look_around(tok,OneLiner);
            Objective=[Objective,[look_behind,tok]];
            if isempty(tok)
                Objective=[Objective,OneLiner];
            end
        end
        OneLiner=rest;
    end
    if strcmp(Objective(1),';')
        error([mfilename,':: No objective function provided. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
    end
    tmp={line_number_,Objective,fichier_name_; % policy objective
        line_number_,Commitment,fichier_name_; % commitment
        line_number_,Discount,fichier_name_};  % discount
    for ii=1:size(tmp,1)
        PlannerObjective_block=capture_equations(PlannerObjective_block,tmp(ii,:),'planner_objective');
    end
    % now make sure there are variables and there is no time in the variables
    time_vars=cell2mat(PlannerObjective_block{1}(2,~cellfun(@isempty,PlannerObjective_block{1}(2,:))));
    if isempty(time_vars)
        error([mfilename,':: Planner objective must include variables. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
    elseif any(time_vars)
        error([mfilename,':: leads or lags not allowed in planner_objective block. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
    end
    dictionary.planner.model={Objective;Commitment;Discount};
end

% remove item from block
blocks(current_block_id)=[];

%% parameterization block
current_block_id=find(strcmp('parameterization',{blocks.name}));
% % if ~is_svar_modelisempty(blocks(current_block_id).listing)
% %     warning([mfilename,':: parameterization block missing:: may not be able to take the derivatives wrt switching parameters...']) %#ok<WNTAG>
% % end
dictionary.Parameterization_block=cell(0,1);
% first make sure all statements are on the same line
blocks(current_block_id).listing=reprocess_parameter_batch(blocks(current_block_id).listing);
% now do it
for ii=1:size(blocks(current_block_id).listing,1)
    dictionary.Parameterization_block=capture_parameterization(dictionary.Parameterization_block,blocks(current_block_id).listing(ii,:));
end
% remove item from block
blocks(current_block_id)=[];

% check that every parameter is controled by one chain only and that it is
% assigned a value in every state of the commanding chain. NB: above, I
% have already made sure that every parameter in the parameterization block
% is declared
for ii=1:numel(dictionary.parameters)
    par_i=dictionary.parameters(ii).name;
    loc=find(strcmp(par_i,dictionary.Parameterization_block(:,1)));
    [istp,diagonal]=is_transition_probability(par_i);
    if istp && diagonal
        error([mfilename,':: remove transition probability ''',par_i,''' and list off-diagonal elements only'])
    end
    if ~isempty(loc)
        chains=unique(dictionary.Parameterization_block(loc,2));
        if numel(chains)>1
            disp(chains(:)')
            error([mfilename,':: parameter ''',par_i,''' controlled by multiple chains'])
        end
        
        % check that there no duplicate statements
        for j1=1:numel(loc)
            for j2=j1+1:numel(loc)
                if isequal(dictionary.Parameterization_block{loc(j1),3},dictionary.Parameterization_block{loc(j2),3})
                    error([mfilename,':: parameterization of ',par_i,' occurs at least twice in the same state'])
                end
            end
        end
        
        % check the each parameter is assigned a value in each state of the
        % governing chain
        states_occur=cell2mat(dictionary.Parameterization_block(loc,3));
        chain_loc=strcmp(chains,dictionary.MarkovChains(1,:));
        chain_states=1:dictionary.MarkovChains{2,chain_loc};
        missing=~ismember(chain_states,states_occur);
        if any(missing)
            error([mfilename,':: parameter ''',par_i,...
                ''' not assigned a value in states ',...
                num2str(chain_states(missing))])
        end
        missing=~ismember(states_occur,chain_states);
        if any(missing)
            error([mfilename,':: parameter ''',par_i,...
                ''' assigns values in states ',...
                num2str(states_occur(missing)),' but those states are not declared in chain ',chains])
        end
        if numel(states_occur)>1
            dictionary.parameters(ii).is_switching=true;
        end
    end
end

% check that all the chains are fully described. This does not concern
% constant parameters
for ii=2:size(dictionary.MarkovChains,2)
    chain_name=dictionary.MarkovChains{1,ii};
    if dictionary.MarkovChains{3,ii} % endogenous
        membership=dictionary.time_varying_probabilities;
    else
        membership={dictionary.parameters.name};
    end
    if dictionary.MarkovChains{2,ii}<2 % number of states
        error([mfilename,':: Markov chain ''',chain_name,''' must have at least 2 states. The one state case is for constant parameters only'])
    else
        missing=cell(1,0);
        for jj=1:dictionary.MarkovChains{2,ii}
            for kk=1:dictionary.MarkovChains{2,ii}
                prob_jk=[chain_name,'_tp_',int2str(jj),'_',int2str(kk)];
                if kk~=jj && ~ismember(prob_jk,membership)
                    missing=[missing,{prob_jk}];
                end
            end
        end
        if ~isempty(missing)
            disp(missing)
            error([mfilename,':: the chains just listed are missing and must be declared'])
        end
    end
end

%% parameter restrictions block.

current_block_id=find(strcmp('parameter_restrictions',{blocks.name}));
Param_rest_block=cell(0,1);
for ii=1:size(blocks(current_block_id).listing,1)
    Param_rest_block=capture_equations(Param_rest_block,blocks(current_block_id).listing(ii,:),'parameter_restrictions');
end
% remove item from block
blocks(current_block_id)=[]; %#ok<NASGU>

dictionary.Param_rest_block=Param_rest_block;
clear Param_rest_block

%% Lump together the model and steady-state model
AllModels=[Model_block;
    SteadyStateModel_block;
    PlannerObjective_block];
ss_eq_nbr=numel(SteadyStateModel_block);
% steady state equations are identified by number 4
planobj_eq_nbr=numel(PlannerObjective_block);
% dictionary.planner objective equations are identified by number 5
equation_type=[equation_type;4*ones(ss_eq_nbr,1);5*ones(planobj_eq_nbr,1)];

clear Model_block SteadyStateModel_block PlannerObjective_block
%% replace log variables and create new variables
% old endo names will be useful if there is a problem in the model and has
% to be constructed prior to changing the names of the variables.
old_endo_names={dictionary.orig_endogenous.name};
if numel(dictionary.log_vars)
    logvarnames={dictionary.log_vars.name};
    % check the model is not declared to be linear
    if dictionary.is_linear_model
        error([mfilename,':: with ''log_vars'', the model cannot be declared as ''linear'''])
    end
    % re-sweep all the equations
    for ii=1:numel(AllModels)
        eq_i=AllModels{ii};
        new_eq_i=cell(2,0);
        for jj=1:size(eq_i,2)
            new_item=eq_i(:,jj);
            if ~isempty(eq_i{2,jj}) && ismember(eq_i{1,jj},logvarnames)
                new_item{1,1}=['LOG_',new_item{1,1}];
                % check whether we are dealing with a steady state first
                if size(new_eq_i,2)>=2 && strcmp(new_eq_i{1,end-1},'steady_state')
                    new_eq_i=new_eq_i(:,1:end-2);
                    % open two parentheses but close only one
                    new_item=[{'exp(',[]}',{'steady_state',[]}',{'(',[]}',new_item,{')',[]}'];
                else
                    new_item=[{'exp(',[]}',new_item,{')',[]}'];
                end
            end
            new_eq_i=[new_eq_i,new_item];
        end
        % before proceeding, check whether it is an equation with equality,
        % in which case we might be dealing with the steady state.
        if size(new_eq_i,2)>4 && strcmp(new_eq_i{1,1},'exp(') && strcmp(new_eq_i{1,4},'=')
            % locate the semicolon
            if strcmp(new_eq_i{1,end},';')
                new_eq_i=new_eq_i(:,1:end-1);
            else
                new_eq_i{1,end}=new_eq_i{1,end}(1:end-1);
            end
            % invert the equation and add the suppressed semicolon
            new_eq_i=[new_eq_i(:,2),new_eq_i(:,4),{'log(',[]}',new_eq_i(:,5:end),{');',[]}'];
        end
        AllModels{ii}=new_eq_i;
    end
    % then change the names in the list of endogenous
    for ii=1:numel(dictionary.log_vars)
        vname=dictionary.log_vars(ii).name;
        var_pos=strcmp(vname,{dictionary.orig_endogenous.name});
        dictionary.orig_endogenous(var_pos).name=['LOG_',vname];
    end
    % re-order the variables
    [~,tags]=sort({dictionary.orig_endogenous.name});
    dictionary.orig_endogenous=dictionary.orig_endogenous(tags);
    % re-order the incidence and occurrence matrices accordingly
    Occurrence=Occurrence(:,tags,:);
end

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
% % % % % % % % % % dictionary.Lead_lag_incidence=transpose(dictionary.Lead_lag_incidence);

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
dictionary.planner.shadow_model=cell(0,1);

static.steady_state_shadow_model=cell(0,1);
shadow_definitions=cell(0,1);
orig_definitions=cell(0,1);
shadow_tvp=cell(0,1);

for ii=1:numel(equation_type)
    % we don't need the semicolons anymore
    eq_i=AllModels{ii};
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
        [status,pos]=determine_status(item);
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
                    sh_tvp=[sh_tvp,'y(',int2str(index),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'y(',int2str(index),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'y(',int2str(index),')'];
                else
                    o_m=[o_m,item];
                    if lead_or_lag
                        o_m=[o_m,'{',int2str(lead_or_lag),'}'];
                    end
                    sh_o=[sh_o,'y(',int2str(index),')'];
                    if ss_is_on
                        sh_vo=[sh_vo,'y(',int2str(index),')'];
                    else
                        sh_vo=[sh_vo,'y(',int2str(index),',:)'];
                    end
                    s_m=[s_m,item];
                    sh_s=[sh_s,'y(',int2str(pos),')'];
                    switch lead_or_lag
                        case -1
                            sh_b1=[sh_b1,'(','y(',int2str(pos),')-y(',int2str(pos+orig_endo_nbr),'))'];
                            sh_b2=[sh_b2,'y(',int2str(pos),')'];
                        case 0
                            sh_b1=[sh_b1,'y(',int2str(pos),')'];
                            sh_b2=[sh_b2,'(','y(',int2str(pos),')+y(',int2str(pos+orig_endo_nbr),'))'];
                        case 1
                            sh_b1=[sh_b1,'(','y(',int2str(pos),')+y(',int2str(pos+orig_endo_nbr),'))'];
                            sh_b2=[sh_b2,'(','y(',int2str(pos),')+2*y(',int2str(pos+orig_endo_nbr),'))'];
                    end
                end
                ss_is_on=false;
            case 'x'
                if is_def
                    error([mfilename,':: definitions cannot contain variables']);
                elseif is_tvp
                    error([mfilename,':: endogenous switching probabilities cannot contain shocks. Use auxiliary variables if necessary']);
                elseif is_sseq
                    sh_ssm=[sh_ssm,'x(',int2str(pos),')'];
                elseif is_planner
                    error([mfilename,':: dictionary.planner objectives cannot include shocks. Use auxiliary variables if necessary'])
                else
                    o_m=[o_m,item];
                    sh_o=[sh_o,'x(',int2str(pos),')'];
                    sh_vo=[sh_vo,'x(',int2str(pos),',:)'];
                    s_m=[s_m,item];
                    sh_s=[sh_s,'x(',int2str(pos),')'];
                    sh_b1=[sh_b1,'x(',int2str(pos),')'];
                    sh_b2=[sh_b2,'x(',int2str(pos),')'];
                end
            case 'param'
                if is_def
                    sh_d=[sh_d,'param(',int2str(pos),')'];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,'param(',int2str(pos),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'param(',int2str(pos),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'param(',int2str(pos),')'];
                else
                    o_m=[o_m,item];
                    sh_o=[sh_o,'param(',int2str(pos),')'];
                    sh_vo=[sh_vo,'param(',int2str(pos),')'];
                    s_m=[s_m,item];
                    sh_s=[sh_s,'param(',int2str(pos),')'];
                    sh_b1=[sh_b1,'param(',int2str(pos),')'];
                    sh_b2=[sh_b2,'param(',int2str(pos),')'];
                end
            case 'def'
                if is_def
                    sh_d=[sh_d,'def(',int2str(pos),')'];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,'def(',int2str(pos),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'def(',int2str(pos),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'def(',int2str(pos),')'];
                else
                    o_m=[o_m,item];
                    s_m=[s_m,item];
                    sh_o=[sh_o,'def(',int2str(pos),')'];
                    sh_vo=[sh_vo,'def(',int2str(pos),')'];
                    sh_s=[sh_s,'def(',int2str(pos),')'];
                    sh_b1=[sh_b1,'def(',int2str(pos),')'];
                    sh_b2=[sh_b2,'def(',int2str(pos),')'];
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
                    eq_i{1,jj+2}=int2str(pos);
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
        dictionary.planner.shadow_model=[dictionary.planner.shadow_model;{sh_pl}];
    else
        dynamic.model=[dynamic.model;{o_m}];
        static.model=[static.model;{s_m}];
        dynamic.shadow_model=[dynamic.shadow_model;{replace_steady_state_call(sh_o)}];
        static.shadow_model=[static.shadow_model;{replace_steady_state_call(sh_s)}];
        static.shadow_BGP_model=[static.shadow_BGP_model;{replace_steady_state_call(sh_b1)};{replace_steady_state_call(sh_b2)}];
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
    static.steady_state_shadow_model=[{['y=zeros(',int2str(orig_endo_nbr),',1);']};static.steady_state_shadow_model];
end

%% symbolic forms

if ~is_svar_model
    % dynamic model wrt y and x
    wrt={'y',1:nnz(dictionary.Lead_lag_incidence)
        'x',1:numel(dictionary.exogenous)};
    [model_derivatives,numEqtns,numVars,jac_toc]=...
        rise_derivatives(dynamic.shadow_model,dictionary.input_list,wrt);
    disp([mfilename,':: first-order derivatives of dynamic model wrt y(+0-) and x. ',...
        int2str(numEqtns),' equations and ',int2str(numVars),' variables :',int2str(jac_toc),' seconds'])
    
    % static model wrt y
    wrt={'y',1:size(dictionary.Lead_lag_incidence,1)};
    [static_model_derivatives,numEqtns,numVars,jac_toc]=...
        rise_derivatives(static.shadow_model,dictionary.input_list,wrt);
    disp([mfilename,':: first-order derivatives of static model wrt y(0). ',...
        int2str(numEqtns),' equations and ',int2str(numVars),' variables :',int2str(jac_toc),' seconds'])
    
    % balanced growth path model wrt y
    wrt={'y',1:2*size(dictionary.Lead_lag_incidence,1)};
    [static_bgp_model_derivatives,numEqtns,numVars,jac_toc]=...
        rise_derivatives(static.shadow_BGP_model,dictionary.input_list,wrt);
    disp([mfilename,':: first-order derivatives of static BGP model wrt y(0). ',...
        int2str(numEqtns),' equations and ',int2str(numVars),' variables :',int2str(jac_toc),' seconds'])
    
    % dynamic model wrt param: I probably should insert the definitions but
    % they are taken out for the moment...
    wrt={'param',1:numel(dictionary.parameters)};
    if DefaultOptions.definitions_in_param_differentiation
        [param_derivatives,numEqtns,numVars,jac_toc]=...
            rise_derivatives(dynamic.shadow_model,dictionary.input_list,wrt,...
            shadow_definitions);
    else
        [param_derivatives,numEqtns,numVars,jac_toc]=...
            rise_derivatives(dynamic.shadow_model,dictionary.input_list,wrt);
    end
    disp([mfilename,':: first-order derivatives of dynamic model wrt param. ',...
        int2str(numEqtns),' equations and ',int2str(numVars),' variables :',int2str(jac_toc),' seconds'])
    
    %% push derivatives in to functions
    % for the moment, put all in the same matrix
    dynamic.model_derivatives=struct(...
        'Endogenous_Shocks',{model_derivatives},...
        'Parameters',{param_derivatives},...
        'StaticEndogenous',{static_model_derivatives},...
        'Static_BGP_Endogenous',{static_bgp_model_derivatives}...
        );
    
    if is_model_with_planner_objective
        wrt={'y',1:size(dictionary.Lead_lag_incidence,1)};
        [LossComDiscHessJac,numEqtns,numVars,jac_toc]=...
            rise_derivatives(dictionary.planner.shadow_model(1),dictionary.input_list,wrt,shadow_definitions,2);
        disp([mfilename,':: 1st and 2nd-order derivatives of planner objective wrt y(0). ',...
            int2str(numEqtns),' equations and ',int2str(numVars),' variables :',int2str(jac_toc),' seconds'])
        % add the loss, the commitment degree and discount
        shadow_model=dictionary.planner.shadow_model;
        lcd={
            ['loss=',shadow_model{1}]
            ['commitment=',strrep(shadow_model{2},'commitment-','')]
            ['discount=',strrep(shadow_model{3},'discount-','')]
            };
        LossComDiscHessJac.code=[cell2mat(lcd(:)'),LossComDiscHessJac.code];
        LossComDiscHessJac.argouts={'loss','commitment','discount','Hess_','Jac_'};
        dictionary.planner.LossComDiscHessJac=LossComDiscHessJac;
    end
end
%% give greek names to endogenous, exogenous, parameters
dictionary.orig_endogenous=greekify(dictionary.orig_endogenous);
dictionary.exogenous=greekify(dictionary.exogenous);
dictionary.parameters=greekify(dictionary.parameters);

%% Add final variables list to the dictionary as they may differ earlier ones
unsorted_endogenous=dictionary.orig_endogenous;

dictionary.NumberOfEquations=sum(equation_type==1);

% finally check that the number of equations is consistent with the number
% of variables

dictionary.is_svar_model=is_svar_model;

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
        new_var=struct('name',['mult_',int2str(eq)],'tex_name','','max_lead',0,'max_lag',0);
        unsorted_endogenous=[unsorted_endogenous,new_var];
    end
elseif ~is_svar_model
    assert(numel(dictionary.orig_endogenous)==sum(equation_type==1),...
        '# equations different from # endogenous variables')
    if dictionary.is_sticky_information_model
        dictionary.forward_looking_ids=find(dictionary.Lead_lag_incidence(:,3));
        for ii=1:numel(dictionary.forward_looking_ids)
            id=dictionary.forward_looking_ids(ii);
            new_var=struct('name',['SI_',dictionary.orig_endogenous{id}],'tex_name','','max_lead',0,'max_lag',0);
            unsorted_endogenous=[unsorted_endogenous,new_var];
        end
    end
end
% now we can resort the final variables
[~,tags]=sort({unsorted_endogenous.name});
dictionary.endogenous=unsorted_endogenous(tags);
% this will be used to reorder the solution of loose commitment or sticky
% information
dictionary.reordering_index=locate_variables({dictionary.endogenous.name},{unsorted_endogenous.name});
dictionary.reordering_index=transpose(dictionary.reordering_index(:));
clear unsorted_endogenous


% greekify the final endogenous
dictionary.endogenous=greekify(dictionary.endogenous);

dictionary.is_endogenous_switching_model=sum(equation_type==3)>0;
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

% locate the multipliers in the ordered list. THIS PROPERTY COULD BE
% IMPLEMENTED DIRECTLY AS A PROPERTY OF THE VARIABLE... BUT IT MIGHT GET
% TOO MESSY EVENTUALLY IF WE STRUCTURE TOO MUCH. I MAY STILL DO THAT IN THE
% FUTURE SINCE THIS PROPERTY IS USED SOMEHOW. IT DEPENDS HOW IT IS USED
dictionary.is_lagrange_multiplier=strncmp('mult_',{dictionary.endogenous.name},5);
dictionary.is_endogenous_switching_model=any([dictionary.MarkovChains{3,:}]);

% load transition matrices
transition_probabilities();

% variable names for the original endogenous but without the augmentation
% add-ons
%dictionary.orig_endo_names_current={orig_endogenous_current.name};
dictionary.dynamic=dynamic;
dictionary.static=static;
clear static dynamic

for ii=1:numel(dictionary.parameters)
    if is_transition_probability(dictionary.parameters(ii).name)
        is_in_use_parameter(ii)=true;
    end
end
dictionary.is_in_use_parameter=is_in_use_parameter;
dictionary.is_in_use_shock=is_in_use_shock;
% replace the list of definition names with definition equations
dictionary.definitions=struct('model',orig_definitions,'shadow',shadow_definitions);
dictionary=orderfields(dictionary);
%% Miscellaneous
% 1.

%% functions

    function transition_probabilities()
        
        ProbScript=[shadow_tvp;{'Q=1;'}];
        shadow_tvp_left_right=cell(1,2);
        for iprob=1:numel(shadow_tvp)
            pp=shadow_tvp{iprob};
            eq_loc=strfind(pp,'=');
            shadow_tvp_left_right{iprob,1}=pp(1:eq_loc-1);
            shadow_tvp_left_right{iprob,2}=pp(eq_loc+1:end-1);
        end
        iter=numel(ProbScript);
        Qsize=1;
        QQ=struct([]); % this will not be used anymore
        for i1=1:size(dictionary.MarkovChains,2)
            chain=dictionary.MarkovChains{1,i1};
            states_nbr_ii=dictionary.MarkovChains{2,i1};
            Qsize=Qsize*states_nbr_ii;
            exo_flag=~logical(dictionary.MarkovChains{3,i1});
            iter=iter+1;
            ProbScript{iter,1}=['Qi=zeros(',int2str(states_nbr_ii),');'];
            NewQ=cell(states_nbr_ii);
            for s1=1:states_nbr_ii
                cumul='0';
                for s2=1:states_nbr_ii
                    if s1~=s2
                        prob_name=[chain,'_tp_',int2str(s1),'_',int2str(s2)];
                        if exo_flag
                            loc=find(strcmp(prob_name,{dictionary.parameters.name}));
                            if isempty(loc)
                                error([mfilename,':: transition probability ',prob_name,' uncharacterized'])
                            end
                            new_tp=['param(',int2str(loc),')'];
                            Info=['Qi(',int2str(s1),',',int2str(s2),')=',new_tp,';'];
                        else
                            prob_loc= strcmp(shadow_tvp_left_right(:,1),prob_name);
                            new_tp=shadow_tvp_left_right{prob_loc,2};
                            Info=['Qi(',int2str(s1),',',int2str(s2),')=',new_tp,';'];
                        end
                        NewQ{s1,s2}=new_tp;
                        cumul=strcat(cumul,'+',new_tp);
                        iter=iter+1;
                        ProbScript{iter,1}=Info;
                    end
                end
                NewQ{s1,s1}=['1-(',cumul,')'];
                iter=iter+1;
                ProbScript{iter,1}=['Qi(',int2str(s1),',',int2str(s1),')=1-sum(Qi(',int2str(s1),',:));'];
            end
            % those probabilities could be functions of
            % definitions. so take care of that right here. 
            NewQ=replace_definitions(NewQ,shadow_definitions);
            QQ(i1).Q=cellstr2mat(NewQ);
            iter=iter+1;
            ProbScript{iter,1}='Q=kron(Q,Qi);';
        end
        dictionary.shadow_transition_matrix=struct('code',cell2mat(ProbScript(:)'),...
            'argins',{dictionary.input_list},... 
            'argouts',{{'Q'}});
%         dictionary.shadow_transition_matrix=QQ;
    end

    function string=replace_steady_state_call(string,flag)
        if nargin<2
            flag='';
        end
        loc_=strfind(string,'steady_state');
        span=length('steady_state');
        while ~isempty(loc_)
            loc_=loc_(1);
            left=string(1:loc_-1);
            right=string(loc_+span:end);
            closing=strfind(right,')');
            closing=closing(1);
            number=right(2:closing-1);
            right=right(closing+1:end);
            switch flag
                case 'symbolic'
                    in_between=['ss_',number];
                otherwise
                    in_between=['ss(',number,')'];
            end
            string=[left,in_between,right];
            loc_=strfind(string,'steady_state');
        end
    end

    function check_validity(syntax,file_name_,iline_,block_name)
        good=ismember(syntax,dictionary.syntax_typical)||...
            (ismember(syntax,dictionary.syntax_time) && ~strcmp(block_name,'steady_state_model'))||...
            (ismember(syntax,dictionary.syntax_function) && function_on)||...
            (ismember(syntax,dictionary.syntax_special) && strcmp(block_name,'steady_state_model'))||...
            (ismember(syntax,dictionary.syntax_special) && strcmp(block_name,'parameter_restrictions'));
        if ~good
            error([mfilename,':: wrong syntax ',syntax,' in ',file_name_,' at line ',int2str(iline_)])
        end
    end

    function [status,loc_]=determine_status(tokk)
        loc_=[];
        possibilities={'orig_endogenous','y'
            'exogenous','x'
            'parameters','param' % same as definitions
            'definitions','def'
            'time_varying_probabilities','tvp'
            'known_words','f' % same as functions
            'symbols',tokk
            'add_operators','+'
            'mult_operators','*'
            'relational_operators','>'
            'chain_names','cn'
            };
        iter=0;
        status='unknown';
        while isempty(loc_) && iter<size(possibilities,1)
            iter=iter+1;
            if iter<4
                loc_=find(strcmp(tokk,{dictionary.(possibilities{iter,1}).name}),1);
            else
                loc_=find(strcmp(tokk,dictionary.(possibilities{iter,1})),1);
            end
            if ~isempty(loc_)
                status=possibilities{iter,2};
                break
            end
        end
        if isempty(loc_)
            if exist([tokk,'.m'],'file')
                status='f'; % same as known words
            else
                try %#ok<TRYNC>
                    isnumeric(eval(tokk));
                    status='n';
                end
            end
        end
        % maybe this function should be called with a block name...
        % %         if strcmp(status,'unknown') && strcmp(current_block_name,'steady_state_model')
        % %             dictionary.known_words=[dictionary.known_words,{tokk}];
        % %             status='f';
        % %         end
    end

    function block=capture_equations(block,cell_info,block_name)
        iline_=cell_info{1};
        rawline_=cell_info{2};
        file_name_=cell_info{3};
        while ~isempty(rawline_) && ~all(isspace(rawline_))
            % try splitting first
            semicol=strfind(rawline_,';');
            if ~isempty(semicol)
                % I don't call look_around here because I need the semicol
                rest_=rawline_(1:semicol(1));
                rawline_=rawline_(semicol(1)+1:end);
            else
                rest_=rawline_;
                rawline_=[];
            end
            
            if isempty(equation)
                % % % % %                 last_status='';
                endo_switch_flag=false;
                def_flag=false;
                if time_on
                    error([mfilename,':: new equation starting without finished time index in ',file_name_,' at line ',int2str(iline_)])
                end
                if function_on
                    error([mfilename,':: new equation starting with earlier function not closed in ',file_name_,' at line ',int2str(iline_)])
                end
            end
            while ~isempty(rest_) && ~all(isspace(rest_))
                [tokk,rest1]=strtok(rest_,DELIMITERS);
                
                % test whether there is a declaration
                if ~isempty(tokk)
                    if (strcmp(block_name,'model') && strcmp(tokk,'model'))||...
                            (strcmp(block_name,'steady_state_model') && strcmp(tokk,'steady_state_model'))||...
                            (strcmp(block_name,'exogenous_definition') && strcmp(tokk,'exogenous_definition'))
                        while ~isempty(rest1)
                            [tokk,rest1]=strtok(rest1,DELIMITERS);
                            if ~isempty(tokk)
                                if strcmp(tokk,'linear') && strcmp(block_name,'model')
                                    dictionary.is_linear_model=true;
                                elseif strcmp(tokk,'imposed') && strcmp(block_name,'steady_state_model')
                                    static.is_imposed_steady_state=true;
                                elseif strcmp(tokk,'unique') && strcmp(block_name,'steady_state_model')
                                    static.is_unique_steady_state=true;
                                else
                                    error([mfilename,':: unknown attribute ''',tokk,''' in file ',file_name_,' at line ',int2str(iline_)])
                                end
                            end
                        end
                        break % exit while ~isempty(rest_) loop
                    end
                end
                
                if ~isempty(tokk)
                    tok_status=determine_status(tokk);
                    if strcmp(tok_status,'param')
                        if strcmp(block_name,'exogenous_definition')
                            error([mfilename,':: exogenous definitions cannot contain parameters in file ',file_name_,' at line ',int2str(iline_)])
                        end
                        position=strcmp(tokk,{dictionary.parameters.name});
                        is_in_use_parameter(position)=true;
                    elseif strcmp(tok_status,'x')
                        position=strcmp(tokk,{dictionary.exogenous.name});
                        is_in_use_shock(position)=true;
                    end
                    if strcmp(tok_status,'#')||strcmp(tok_status,'!')
                        def_flag=strcmp(tok_status,'#');
                        if strcmp(block_name,'exogenous_definition')
                            error([mfilename,':: the exogenous definition block cannot contain ''#'' or ''!'' ',file_name_,' at line ',int2str(iline_)])
                        end
                        endo_switch_flag=strcmp(tok_status,'!');
                        rest_=rest1;
                        [tokk,rest1]=strtok(rest_,DELIMITERS); %#ok<*STTOK>
                        tok_status=determine_status(tokk);
                        if ~strcmp(tok_status,'unknown')
                            if strcmp(tok_status,'f')
                                disp([mfilename,':: (gentle warning): ',tokk,' is also a matlab function'])
                            else
                                error([mfilename,':: string ''',tokk,''' in ',file_name_,' at line ',int2str(iline_),' cannot have multiple types'])
                            end
                        end
                        if def_flag
                            dictionary.definitions=[dictionary.definitions;{tokk}];
                        elseif endo_switch_flag
                            if ~is_transition_probability(tokk)
                                error([mfilename,':: string ''',tokk,''' in ',file_name_,' at line ',int2str(iline_),' is not an appropriate name for an endogenous switching probability'])
                            end
                            dictionary.time_varying_probabilities=[dictionary.time_varying_probabilities,{tokk}];
                        else
                            error([mfilename,':: parsing error in ',file_name_,' at line ',int2str(iline_),' please report this to junior.maih@gmail.com'])
                        end
                        % update the status above
                        tok_status=determine_status(tokk);
                        if ~isempty(equation)
                            error([mfilename,':: # and ! can only occur at the beginning of an equation check in ',file_name_,' line ',int2str(iline_)])
                        end
                    elseif strcmp(tok_status,'unknown')
                        error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',int2str(iline_)])
                    end
                    
                    if def_flag && (strcmp(tok_status,'y')||strcmp(tok_status,'x'))
                        error([mfilename,':: definitions cannot contain variables. check ',file_name_,' at line ',int2str(iline_)])
                    end
                    
                    if ~endo_switch_flag && strcmp(tok_status,'tvp')
                        error([mfilename,':: model equations cannot contain endogenous switching probabilities. check ',file_name_,' at line ',int2str(iline_)])
                    end
                    left_operator=look_around(tokk,rest_);
                    
                    for i1=2:length(left_operator)
                        first=determine_status(left_operator(i1-1));
                        if strcmp(first,'unknown')
                            error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',int2str(iline_)])
                        end
                        second=determine_status(left_operator(i1));
                        if strcmp(second,'unknown')
                            error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',int2str(iline_)])
                        end
                        check_validity([first,second],file_name_,iline_,block_name);
                    end
                    if ~isempty(left_operator)
                        if ~isempty(tokk)
                            third=determine_status(left_operator(end));
                            if strcmp(third,'unknown')
                                error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',int2str(iline_)])
                            end
                            check_validity([third,tok_status],file_name_,iline_,block_name)
                        end
                        % always deal with time before dealing with
                        % functions as time potentially removes some
                        % parentheses and the functions do not.
                        if ((strcmp(last_status,'y')||strcmp(last_status,'x')) &&...
                                (strcmp(left_operator(1),'(')||strcmp(left_operator(1),'{')))% check the time
                            time_on=true;
                            time_opening=left_operator(1);
                            fill_time=[fill_time,left_operator(2:end)];
                            left_operator='';
                        elseif time_on && (strcmp(left_operator(1),')')||strcmp(left_operator(1),'}'))
                            fill_time=eval(fill_time);
                            if ~isnumeric(fill_time)||~isequal(fill_time,floor(fill_time))
                                error([mfilename,':: time syntax error in ',file_name_,' at line ',int2str(iline_)])
                            end
                            equation{2,end}=fill_time;
                            update_leads_lags(equation{1,end},fill_time);
                            fill_time='';
                            time_on=false;
                            time_closing=left_operator(1);
                            if ~((strcmp(time_opening,'(') && strcmp(time_closing,')'))||...
                                    (strcmp(time_opening,'{') && strcmp(time_closing,'}')))
                                error([mfilename,':: time opened with ''',time_opening,''' and closed with ''',time_closing,''' in ',file_name_,' at line ',int2str(iline_)])
                            end
                            left_operator=left_operator(2:end);
                        end
                        
                        left_parents=sum(left_operator=='(');
                        right_parents=sum(left_operator==')');
                        % try and close the function before worrying about
                        % whether the closing is not successful.
                        if function_on && right_parents
                            function_on=function_on-min(function_on,right_parents);
                        end
                        if (strcmp(last_status,'f')||function_on) && left_parents
                            function_on=function_on+left_parents;
                        end
                        if function_on<0
                            error([mfilename,':: parenthesis mismatch in ',file_name_,' at line ',int2str(iline_)])
                        end
                        
                        if ~isempty(left_operator)
                            equation=[equation,{left_operator,[]}'];
                        end
                    end
                    if time_on && strcmp(tok_status,'n')
                        fill_time=[fill_time,tokk];
                    else
                        equation=[equation,{tokk,[]}'];
                        if (strcmp(tok_status,'y')||strcmp(tok_status,'x'))
                            equation{2,end}=0;
                        end
                        fill_time='';
                    end
                    last_status=tok_status;
                    rest_=rest1;
                else
                    if time_on && (strcmp(rest_(1),')')||strcmp(rest_(1),'}'))
                        rest_=rest_(2:end);
                        %====================
                        fill_time=eval(fill_time);
                        if ~isnumeric(fill_time)||~isequal(fill_time,floor(fill_time))
                            error([mfilename,':: time syntax error in ',file_name_,' at line ',int2str(iline_)])
                        end
                        equation{2,end}=fill_time;
                        update_leads_lags(equation{1,end},fill_time);
                        fill_time='';
                        %====================
                        time_on=false;
                    end
                    if function_on
                        left_parents=sum(rest_=='(');
                        right_parents=sum(rest_==')');
                        function_on=function_on+left_parents-min(function_on,right_parents);
                    end
                    if ~isempty(rest_)
                        equation=[equation,{rest_,[]}'];
                        rest_='';
                    end
                    last_status=determine_status(equation{1,end}(end));
                end
            end
            if ~isempty(equation)
                if strcmp(equation{1,end}(end),';')
                    % we've reach the end of the equation, validate it,
                    % load it and reinitialize.
                    equation=validate_equation(equation);
                    block=[block;{equation}];
                    equation=cell(2,0);
                end
            end
        end
        function equation=validate_equation(equation)
            checked=false;
            % look for equality signs
            equality_signs=cell(2,0);
            % the top cell with register the cell number with equality
            % sign and the lower cell, the location of the equality sign
            % within the top cell
            for ic=1:size(equation,2)
                eq_s=strfind(equation{1,ic},'=');
                if ~isempty(eq_s)
                    equality_signs=[equality_signs,transpose({ic,eq_s})];
                end
                if strcmp(block_name,'parameter_restrictions')
                    % the following restrictions must be fulfilled for the
                    % parameter restriction block
                    % 1- no variable
                    if ismember(equation{1,ic},{dictionary.orig_endogenous.name})
                        error([mfilename,':: no variable allowed in the parameter_restrictions block'])
                    end
                    % 2- if parenthesis after a parameter, then one of the elements
                    % inside the parenthesis must be a chain name and the other a
                    % numeric.
                    if ismember(equation{1,ic},{dictionary.parameters.name})
                        % locate the parameter is the parameterization
                        % block and collect its max state and chain
                        % name
                        locs=find(strcmp(equation{1,ic},dictionary.Parameterization_block(:,1)));
                        if isempty(locs)
                            error([mfilename,':: parameter ',equation{1,ic},' not parameterized in ',file_name_,' at line ',int2str(iline_)])
                        end
                        nom_de_chaine=dictionary.Parameterization_block{locs(1),2};
                        etats_de_la_chaine=cell2mat(dictionary.Parameterization_block(locs,3));
                        % put a pointer in lieu of the parameter and its
                        % specification
                        state_=1;
                        chain_='const';
                        if  ic<size(equation,2) && strcmp(equation{1,ic+1}(1),'(')
                            state_=[];
                            chain_=[];
                            if ismember(equation{1,ic+2},dictionary.chain_names) % chain or state
                                if isempty(chain_)
                                    chain_=equation{1,ic+2};
                                else
                                    error([mfilename,':: chain declared more than once for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                                end
                            elseif ~isnan(eval(equation{1,ic+2}))
                                if isempty(state_)
                                    state_=eval(equation{1,ic+2});
                                else
                                    error([mfilename,':: chain declared more than once for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                                end
                            else
                                error([mfilename,':: unrecognized character ',equation{1,ic+2},' for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                            end
                            if ~strcmp(equation{1,ic+3},',') % remove comma
                                error([mfilename,':: comma missing between chain and state for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                            end
                            if ismember(equation{1,ic+4},dictionary.chain_names) % chain or state
                                if isempty(chain_)
                                    chain_=equation{1,ic+4};
                                else
                                    error([mfilename,':: chain declared more than once for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                                end
                            elseif ~isnan(eval(equation{1,ic+4}))
                                if isempty(state_)
                                    state_=eval(equation{1,ic+4});
                                else
                                    error([mfilename,':: chain declared more than once for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                                end
                            end
                            if ~strcmp(equation{1,ic+5}(1),')') % closing parenthesis
                                error([mfilename,':: closing parenthesis missing for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                            end
                            % remove all processed elements
                            equation{1,ic+1}(1)=[]; % left parenthesis
                            equation{1,ic+2}=''; % state or chain
                            equation{1,ic+3}=''; % comma
                            equation{1,ic+4}=''; % state or chain
                            equation{1,ic+5}(1)=[]; % right parenthesis
                        end
                        if ~strcmp(chain_,nom_de_chaine)||~(isnumeric(state_) && ismember(state_,etats_de_la_chaine))
                            error([mfilename,':: chain or state missing or wrongly specified for parameter ',equation{1,ic},' in ',file_name_,' at line ',int2str(iline_)])
                        end
                        % we carefully record the chain and the state for
                        % later processing
                        equation{2,ic}={chain_,state_};
                    end
                    % 3- parameter must be effectively controled by the
                    % chain in the parameterization block
                end
            end
            % an equation cannot have multiple equality signs
            if size(equality_signs,2)>1
                error([mfilename,':: multiple equality signs found in ',file_name_,' at line ',int2str(iline_)])
            elseif size(equality_signs,2)==1
                % validate the lhs and the rhs of the equality separately
                lhs=equation(:,1:equality_signs{1}-1);
                rhs=equation(:,equality_signs{1}+1:end);
                % split the middle cell and isolate the equality sign
                string=equation{1,equality_signs{1}};
                left_string=string(1:equality_signs{2}-1);
                right_string=string(equality_signs{2}+1:end);
                if ~isempty(left_string)
                    % add to the lhs
                    lhs=[lhs,transpose({left_string,[]})];
                end
                msg=count_parentheses(lhs(1,:));
                if ~isempty(msg)
                    error([mfilename,':: ',msg,' on left hand side in ',file_name_,' at line ',int2str(iline_)])
                end
                
                if ~isempty(right_string)
                    % add to the rhs but before !
                    rhs=[transpose({right_string,[]}),rhs];
                end
                msg=count_parentheses(rhs(1,:));
                if ~isempty(msg)
                    error([mfilename,':: ',msg,' on right hand side in ',file_name_,' at line ',int2str(iline_)])
                end
                checked=true;
                % Putting the equation back together: what to do with the
                % equality sign?
                if ~strcmp(block_name,'parameter_restrictions') &&...
                        ~strcmp(block_name,'steady_state_model') &&...
                        ~strcmp(block_name,'exogenous_definition') &&...
                        ~ def_flag &&... % <--- should do the same as here with TVP below
                        ~ismember(equation{1,1},dictionary.time_varying_probabilities)
                    % % % % %                         ~ismember(equation{1,1},dictionary.known_words) &&... % I don't remember why this is here !
                    % modify the last item of the right hand side
                    rhs{1,end}=[rhs{1,end}(1:end-1),');'];
                    middle=transpose({'-(',[]});
                else
                    middle=transpose({'=',[]});
                end
                equation=[lhs,middle,rhs];
            end
            if ~checked
                % then there was no equality sign and so, check the whole
                % equation in one go
                msg=count_parentheses(equation(1,:));
                if ~isempty(msg)
                    error([mfilename,':: ',msg,' in ',file_name_,' at line ',int2str(iline_)])
                end
            end
            % if it is an exogenous definition equation, then it needs a
            % special treatment
            if strcmp(block_name,'exogenous_definition')
                % check that the first element on the left-hand side is an observable
                % exogenous
                vnamex=equation{1,1};
                if ~(ismember(vnamex,{dictionary.exogenous.name}) && ...
                        ismember(vnamex,{dictionary.observables.name}))
                    error([mfilename,':: ',vnamex,' must be exogenous and observed in ',file_name_,' at line ',int2str(iline_)'])
                end
                % It must be dated at time 0
                if ~isequal(equation{2,1},0)
                    error([mfilename,':: ',vnamex,' must be dated at time 0 in ',file_name_,' at line ',int2str(iline_)'])
                end
                % the second element must be an equality sign
                if ~strcmp(equation{1,2}(1),'=')
                    error([mfilename,':: the second token must be a ''=''  in ',file_name_,' at line ',int2str(iline_)])
                end
                % no rhs variable can be dated at a time beyond 0
                for icol=3:size(equation,2)
                    if ~isempty(equation{2,icol}) && equation{2,icol}>0
                        error([mfilename,':: right-hand side variables cannot be dated in the future in ',file_name_,' at line ',int2str(iline_)'])
                    end
                end
            end
            function msg=count_parentheses(equation)
                equation=cell2mat(equation);
                leftpars=numel(strfind(equation,'('));
                rightpars=numel(strfind(equation,')'));
                checked_=(leftpars==rightpars);
                leftpars=numel(strfind(equation,'['));
                rightpars=numel(strfind(equation,']'));
                checked_=checked_+2*(leftpars==rightpars);
                switch checked_
                    case 3
                        msg='';
                    case 2
                        msg='parentheses mismatch';
                    case 1
                        msg='brackets mismatch';
                    case 0
                        msg='parentheses and brackets mismatch';
                end
            end
        end
        function update_leads_lags(variable,leadorlag)
            position=find(strcmp(variable,{dictionary.orig_endogenous.name}));
            if ~isempty(position)
                dictionary.orig_endogenous(position).max_lag=min(dictionary.orig_endogenous(position).max_lag,leadorlag);
                dictionary.orig_endogenous(position).max_lead=max(dictionary.orig_endogenous(position).max_lead,leadorlag);
            else
                position=find(strcmp(variable,{dictionary.exogenous.name}));
                if ~isempty(position)
                    if leadorlag>0
                        error([mfilename,':: exogenous variable ',equation{1,end},' cannot have leads. Check ',file_name_,' at line ',int2str(iline_)])
                    end
                    dictionary.exogenous(position).max_lag=min(dictionary.exogenous(position).max_lag,leadorlag);
                end
            end
        end
    end

    function newbatch=reprocess_parameter_batch(batch)
        % this is to make sure all complete statements are on the same line
        newbatch=cell(0,3);
        lines=[];
        dough='';
        while ~isempty(batch)
            bb=batch(1,:);
            lines=[lines,bb{1}];
            dough=[dough,bb{2}];
            dough=strtrim(dough);
            if strcmp(dough(end),';')
                newbatch=[newbatch;
                    {lines,dough,bb{3}}];
                lines=[];
                dough='';
            end
            batch=batch(2:end,:);
        end
    end

    function block=capture_parameterization(block,cell_info)
        PARAM_DELIMITERS=[char([9:13,32]),',;'];
        iline_=cell_info{1};
        old_line=iline_;
        iline_=int2str(old_line(1));
        for iii=2:numel(old_line)
            iline_=[iline_,' & ',int2str(old_line(iii))];
        end
        rawline_=cell_info{2};
        rawline_(isspace(rawline_))=[];
        % remove the semi-colon
        if strcmp(rawline_(end),';')
            rawline_=rawline_(1:end-1);
        end
        file_name_=cell_info{3};
        if isempty(rawline_)
            return
        end
        % name, chain, state, startval, plb, pub, distr      , prob
        % nan,  const,    1 ,    nan  , nan, nan, uniform_pdf,  1
        parameter_array=cell(1,8);
        [tokk,rawline_]=strtok(rawline_,[PARAM_DELIMITERS,'(']);
        if ismember(tokk,{dictionary.observables.name}) && ~ismember(tokk,{dictionary.exogenous.name})
            % create a new parameter and sort the list. At this point,
            % sorting is safe if we have not started rewriting the model in
            % terms of y, x, param, def. Furthermore, there will be or
            % there has already been a check that all observables are
            % either endogenous or exogenous... Only endogenous observables
            % can have a measurement error.
            tokk=['stderr_',tokk];
            newparam=struct('name',tokk,'tex_name','','max_lead',0,'max_lag',0,'is_switching',false,'is_measurement_error',true);
            dictionary.parameters=[dictionary.parameters,newparam];
            clear newparam
            [~,tags]=sort({dictionary.parameters.name});
            dictionary.parameters=dictionary.parameters(tags);
            % the new parameter is of course in use.
            is_in_use_parameter=[is_in_use_parameter;true];
            is_in_use_parameter=is_in_use_parameter(tags);
        end
        if ~ismember(tokk,{dictionary.parameters.name})
            error([mfilename,':: ',tokk,' is not recognized as a parameter or as an observable variable in ',file_name_,' at line(s) ',iline_])
        end
        % diagonal transition probabilities not estimated
        [~,diagonal]=is_transition_probability(tokk);
        if diagonal
            error([mfilename,':: diagonal transition probabilities are not allowed in ',file_name_,' at line(s) ',iline_])
        end
        parameter_array{1}=tokk;
        if strcmp(rawline_(1),'(')
            right_par=strfind(rawline_,')');
            right_par=right_par(1);
            extract=rawline_(2:right_par-1);
            rawline_=rawline_(right_par+1:end);
            if ~isempty(extract)
                % locate the comma and split
                comma=strfind(extract,',');
                if isempty(comma)||numel(comma)~=1
                    error([mfilename,':: exacly 2 arguments are to enter the parentheses right after the parameter name in ',file_name_,' at line(s) ',iline_])
                else
                    first=extract(1:comma-1);
                    second=extract(comma+1:end);
                    if ismember(first,dictionary.MarkovChains(1,:))
                        chain=first;
                        try
                            eval_tok=eval(second);
                        catch %#ok<*CTCH>
                            error([mfilename,':: parameter ',parameter_array{1},' has nonsensical expression for its state in ',file_name_,' at line(s) ',iline_])
                        end
                        state=eval_tok;
                    else
                        try
                            eval_tok=eval(first);
                        catch %#ok<*CTCH>
                            error([mfilename,':: parameter ',parameter_array{1},' has nonsensical expression for its state in ',file_name_,' at line(s) ',iline_])
                        end
                        state=eval_tok;
                        if ismember(second,dictionary.MarkovChains(1,:))
                            chain=second;
                        else
                            error([mfilename,':: ',second,' not in the list of chain names in ',file_name_,' at line(s) ',iline_])
                        end
                    end
                end
            else
                error([mfilename,':: content of parentheses after the parameter name should not be empty in ',file_name_,' at line(s) ',iline_])
            end
        else
            chain='const';
            state=1;
        end
        parameter_array{2}=chain;
        parameter_array{3}=state;
        iter=3;
        
        while ~isempty(rawline_) && iter<8
            [tokk,rawline_]=strtok(rawline_,PARAM_DELIMITERS);
            if ~isempty(tokk)
                iter=iter+1;
                switch iter
                    case {4,5,6,8}
                        try
                            eval_tok=eval(tokk);
                        catch %#ok<*CTCH>
                            error([mfilename,':: lower or upper bound(probability) must be a number. Check ',file_name_,' at line(s) ',iline_])
                        end
                        if ~isreal(eval_tok)
                            error([mfilename,':: lower or upper bound(probability) must be real. Check ',file_name_,' at line(s) ',iline_])
                        end
                        if iter==8 && ~(eval_tok>0||eval_tok<=1)
                            error([mfilename,':: probability must be a real number in (0,1]. Check ',file_name_,' at line(s) ',iline_])
                        end
                    case 7
                        eval_tok=tokk;
                        if ~ismember(eval_tok,dictionary.Distributions)
                            disp(dictionary.Distributions)
                            error([mfilename,':: ',eval_tok,' in ',file_name_,' at line ',int2str(iline_),' is an unknown distribution'])
                        end
                end
                parameter_array{iter}=eval_tok;
            end
        end
        if ~isempty(rawline_)
            error([mfilename,':: apparently more arguments than needed in ',file_name_,' at line(s) ',iline_])
        end
        % now check for errors
        if isempty(parameter_array{4})
            error([mfilename,':: to few arguments for parameterization in ',file_name_,' at line(s) ',iline_])
        end
        
        if isempty(parameter_array{5})
            parameter_array([5,6])={nan};
        end
        
        if isempty(parameter_array{7})
            parameter_array{7}='uniform_pdf';
        end
        
        if isempty(parameter_array{8})
            parameter_array{8}=1;
        end
        
        block=[block;parameter_array];
    end

    function block=construct_list(block,rawline_,tokk)
        if is_trigger
            [~,~,rawline_]=look_around(tokk,rawline_);
        end
        end_game=~isempty(strfind(rawline_,';'));
        while ~isempty(rawline_)
            [tokk,rest_]=strtok(rawline_,DELIMITERS);
            if ~isempty(tokk)
                if strcmp(tokk(1),'$')
                    if ~dollar_active
                        dollar_active=true;
                        tokk=tokk(2:end);
                    end
                end
                if dollar_active
                    second_dollar=strfind(tokk,'$');
                    if ~isempty(second_dollar)
                        tokk=tokk(1:second_dollar-1);
                    end
                    if isempty(tex_name)
                        tex_name=tokk;
                    else
                        if isempty(tokk)
                            if ~isempty(second_dollar)
                                % we've just deleted a dollar sign and so, look
                                % around a $
                                look_behind=look_around('$',rawline_,true);
                            else
                                error([mfilename,':: parsing error. Please contact junior.maih@gmail.com'])
                            end
                        else
                            look_behind=look_around(tokk,rawline_,true);
                        end
                        % preserve space for tex names
                        tex_name=[tex_name,look_behind,tokk]; %#ok<*AGROW>
                        tex_name=strtrim(tex_name);
                    end
                    if ~isempty(second_dollar)
                        dollar_active=false;
                        if ~isempty(block.listing{end,2})
                            error([mfilename,':: found two consecutive tex names in ',...
                                file_name,' at line ',int2str(line_number)])
                        end
                        block.listing{end,2}=tex_name;
                        tex_name='';
                    end
                else
                    if ~isvarname(tokk)
                        error([mfilename,':: atom ''',tokk,''' in ',...
                            file_name,' at line ',int2str(line_number),' is not a valid variable or parameter name'])
                    end
                    if exist([tokk,'.m'],'file')
                        disp([mfilename,':: (gentle warning): ',tokk,' is also a matlab function'])
                    end
                    if ismember(tokk,block.listing(:,1))
                        error([mfilename,':: atom ''',tokk,''' has been declared twice in ',...
                            file_name,' at line ',int2str(line_number)])
                    end
                    block.listing=[block.listing;{tokk,''}];
                end
            end
            rawline_=rest_;
        end
        if end_game
            [~,tags]=sort(block.listing(:,1));
            block.listing=block.listing(tags,:);
            block.active=false;
        end
    end

end
function [is_trigger,blocks,last_block_id]=check_block(tokk,blocks,last_block_id,file_name,line_number)
loc_=find(strcmp(tokk,{blocks.trigger}));
is_trigger=false;
if ~isempty(loc_)
    is_trigger=true;
    if any([blocks.active])
        error([mfilename,':: block ''',blocks(last_block_id).name,...
            ''' active while block ''',blocks(loc_).name ,...
            ''' wants activation in ',file_name,...
            ' at line ',int2str(line_number)])
    end
    blocks(loc_).active=true;
    last_block_id=loc_;
elseif ~any([blocks.active])
    last_block_id=[];
end
end

function [look_behind,loc_,look_forward]=look_around(tokk,to_process,preserve_space)
if nargin<3
    preserve_space=false;
end
loc_=strfind(to_process,tokk);
look_forward=to_process;
look_behind='';
if ~isempty(loc_)
    loc_=loc_(1);
    look_behind=to_process(1:loc_-1);
    if ~preserve_space
        look_behind(isspace(look_behind))=[];
    end
    span=length(tokk);
    look_forward=to_process(loc_+span:end);
end
end


function equation=greekify(equation)
% so far I will just greekify names but later on, I might also greekify
% equations
greek_letters={'alpha','beta','gamma','delta','epsilon','kappa',...
    'lambda','mu','nu','omega','phi','pi','chi','psi','rho',...
    'sigma','tau','upsilon','Sigma','Pi','Lambda','Omega','Gamma'};
DELIMITERS=[char([9:13,32]),'[]{}(),;=+-*/^@'];
% 1-\beta \frac{\left( 1-\frac{\kappa }{2}\left( \Pi _{t}-1\right) ^{2}\right)
% Y_{t}}{\left( 1-\frac{\kappa }{2}\left( \Pi _{t+1}-1\right) ^{2}\right)
% Y_{t+1}}\frac{1}{\exp \left( \mu _{t+1}\right) }\frac{R_{t}}{\Pi _{t+1}}
for ii=1:numel(equation)
    if isempty(equation(ii).tex_name)
        equation(ii).tex_name=greekify_intern(equation(ii).name);
    end
end
    function greek=greekify_intern(equation)
        greek='';
        while ~isempty(equation)
            [tokk,rest_]=strtok(equation,DELIMITERS);
            if isempty(tokk)
                greek=[greek,equation];
            else
                loc_=strfind(equation,tokk);
                loc_=loc_(1);
                greek=[greek,equation(1:loc_-1)];
                span=length(tokk);
                for i1=1:numel(greek_letters)
                    g=greek_letters{i1};
                    gspan=length(g);
                    if span>=gspan
                        if strcmp(tokk(1:span),g) && (span==gspan||strcmp(tokk(span+1),'_'))
                            right=tokk(span+2:end);
                            tokk=['\',g];
                            if ~isempty(right)
                                tokk=[tokk,' _{',right,'}']; % NB: there is a space
                            end
                        end
                    end
                end
                greek=[greek,tokk];
            end
            equation=rest_;
        end
    end
end

% function next=find_next_atom(string,atom)
% next=strfind(string,atom);
% if ~isempty(next)
%     next=next(1);
% else
%     next=0;
% end
% end

% function obj=itemize(obj,nextitem,next_tex_name)
% if nargin<2
%     next_tex_name=[];
% end
% if isempty(obj)
%     obj.list={nextitem};
%     obj.tex_list={next_tex_name};
%     obj.number=1;
% else
%     obj.list=[obj.list;{nextitem}];
%     obj.tex_list=[obj.tex_list;{next_tex_name}];
%     obj.number=obj.number+1;
%     [~,tag]=sort(obj.list);
%     obj.list=obj.list(tag);
%     obj.tex_list=obj.tex_list(tag);
% end
%
% end