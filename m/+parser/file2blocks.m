function [blocks,markov_chains]=file2blocks(RawFile)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

tex_name='';
quote_active=false;
last_block_id=[];

block_col_names={'name','trigger','listing'};

current_markov_chain_name='';
current_number_of_states=[];
new_markov_chain_tex_names={};
markov_chains=parser.initialize_markov_chain('const',1,'is_endogenous',false);

all_chain_names={markov_chains.name};

blocks=parser.initialize_blocks();

% 'endogenous','exogenous','parameters','observables' blocks are just
% declarations. The corresponding blocks are to hold name and tex_name
% respectively

% construct cells with 3 columns that will hold line_number, equation,
% dictionary.filename respectively

blocks=cell2struct(blocks,block_col_names,2);

trigger_map=cell(0,3);
for ii=1:numel(blocks)
    this_triggers=blocks(ii).trigger;
    blknm=blocks(ii).name;
    if ischar(this_triggers)
        this_triggers={this_triggers};
    end
    for itrig=1:numel(this_triggers)
        trigger_map=[trigger_map;{blknm,this_triggers{itrig},ii}];
    end
end

DELIMITERS=parser.delimiters();

blknames={blocks.name};

NumberOfLines=size(RawFile,1);
iline=0;

left_parenth_4_markov_chains_open=false;

while iline<NumberOfLines
    iline=iline+1;
    rawline=RawFile{iline,1};
    file_name=RawFile{iline,2};
    line_number=RawFile{iline,3};
    
    [tok,rest]=strtok(rawline,DELIMITERS);
    if strcmp(tok,'end')
        if isempty(last_block_id)
            error('''end'' found before an opening of a block')
        end
        warning(['ending block ',blocks(last_block_id).name,' with an ''end'' is no longer required. In ',...
            file_name,' at line ',sprintf('%0.0f',line_number)])
        continue
    end
    [is_trigger,last_block_id]=check_block(tok,trigger_map,last_block_id);
    current_block_name=blocks(last_block_id).name;
    switch current_block_name
        case 'Legend'
            % do nothing
        case {'log_vars','endogenous','exogenous','parameters','observables',...
                'level_variables'}
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
            blocks(last_block_id).listing=[blocks(last_block_id).listing;{line_number,rawline,file_name}];
    end
end

% sort the declaration blocks
%----------------------------
sort_list={'log_vars','endogenous','exogenous','parameters','observables',...
    'level_variables'};
for isort=1:numel(sort_list)
    loc=strcmp(sort_list{isort},blknames);
    [~,tags]=sort({blocks(loc).listing.name});
    blocks(loc).listing=blocks(loc).listing(tags);
end

% sort the markov chains
%-----------------------
[~,tags]=sort({markov_chains.name});
markov_chains=markov_chains(tags);

levelVarBlock=strcmp('level_variables',{blocks.name});
levelvar_names={blocks(levelVarBlock).listing.name};
logvarBlock=strcmp('log_vars',{blocks.name});
logvar_names={blocks(logvarBlock).listing.name};
endogBlock=strcmp('endogenous',{blocks.name});
endovar_names={blocks(endogBlock).listing.name};
exogBlock=strcmp('exogenous',{blocks.name});
exovar_names={blocks(exogBlock).listing.name};
obsBlock=strcmp('observables',{blocks.name});
obsVars = {blocks(obsBlock).listing.name};
paramBlock=strcmp('parameters',{blocks.name});
param_names={blocks(paramBlock).listing.name};
modelBlock=strcmp('model',{blocks.name});

if ~isempty(levelvar_names) 
    
    if ~isempty(logvar_names)
        
        error('log_variables and level_variables cannot be declared in the same file')
    
    end
    
    locs=locate_variables(levelvar_names,endovar_names,true);
    
    store_locs=locs;
    
    locs=find(isnan(locs));
    
    if ~isempty(locs)
        
        bad_vars=levelvar_names(locs);
        
        disp(bad_vars(:)')
        
        error('The LEVEL variables above have not been found in the list of endogenous variables')
        
    end
    
    blocks(levelVarBlock).listing=blocks(logvarBlock).listing;
    
%     levelvar_names={blocks(logvarBlock).listing.name};
    
    % now swap and destroy
    blocks(logvarBlock).listing=blocks(endogBlock).listing;
    
    blocks(logvarBlock).listing(store_locs)=[];
    
    logvar_names={blocks(logvarBlock).listing.name};    
    
end

% check that the potential measurement errors have corresponding observables
%--------------------------------------------------------------------------
meas_errs=strncmp('stderr_',param_names,7);
if ~isempty(meas_errs)
    meas_errs=param_names(meas_errs);
    for ii=1:numel(meas_errs)
        obsname=meas_errs{ii}(7+1:end);
        if ~ismember(obsname,obsVars)
            error(['no observable variable with name "',obsname,'" found: "',meas_errs{ii},'" cannot be a measurement error'])
        end
        loc=locate_variables(meas_errs{ii},param_names);
        blocks(paramBlock).listing(loc).is_in_use=true;
        blocks(paramBlock).listing(loc).is_measurement_error=true;
    end
end

% check that all log_vars are endogenous
%---------------------------------------
if ~isempty(logvar_names)
    locs=locate_variables(logvar_names,endovar_names,true);
    locs=find(isnan(locs));
    if ~isempty(locs)
        logvar_names=logvar_names(locs);
        disp(logvar_names(:)')
        error('The LOG variables above have not been found in the list of endogenous variables')
    end
end

% check that all observables are either endogenous or exogenous and are not
% log_vars
%--------------------------------------------------------------------------
for iobs=1:numel(obsVars)
    if ismember(obsVars{iobs},endovar_names)
        blocks(obsBlock).listing(iobs).is_endogenous=true;
    elseif ~ismember(obsVars{iobs},exovar_names)
        error(['observable variable ',obsVars{iobs},' must be either endogenous or exogenous'])
    end
    if ismember(obsVars{iobs},logvar_names)
        error(['observable variable ',obsVars{iobs},' cannot be declared as log_var: use auxiliary variables if necessary'])
    end
end

% check that the same name is not found in other places
%------------------------------------------------------
% no endogenous is exogenous or parameter
for iendo=1:numel(endovar_names)
    if ismember(endovar_names{iendo},exovar_names)
        error(['endogenous variable ',endovar_names{iendo},' cannot be exogenous'])
    elseif ismember(endovar_names{iendo},param_names)
        error(['endogenous variable ',endovar_names{iendo},' cannot be parameter'])
    end
end
% no exogenous is parameter
for iexo=1:numel(exovar_names)
    if ismember(exovar_names{iexo},param_names)
        error(['exogenous variable ',exovar_names{iexo},' cannot be parameter'])
    end
end

% sort the markov chains and mark the parameters, replacing the controlling
% chain with a number representing the order of the markov chain.
% all exogenous transition probabilities are associated with a declared
% markov chain. Later on, it will have to be the case for endogenous
% probabilities as well.
%-----------------------------------------------------------------------
[~,tags]=sort({markov_chains.name});
markov_chains=markov_chains(tags);
markov_chain_names={markov_chains.name};
const_loc=find(strcmp('const',markov_chain_names));
for iparam=1:numel(param_names)
    [istp,~,chain_name,max_state]=parser.is_transition_probability(param_names{iparam});
    if istp
        loc=find(strcmp(chain_name,markov_chain_names));
        if isempty(loc)
            error(['markov chain ',chain_name,' has not been declared'])
        end
        this_nstates=markov_chains(loc).number_of_states;
        if max_state>this_nstates
            error(['maximum number of states of chain "',...
                chain_name,'" was declared to be ',sprintf('%0.0f',...
                this_nstates),', which is inconsistent with ',param_names{iparam}])
        end
        % the transition probabilities are automatically in used
        blocks(paramBlock).listing(iparam).is_in_use=true;
        % so far the chain are all exogenous. This might change when we
        % parse the model equations later on
        blocks(paramBlock).listing(iparam).is_trans_prob=true;
        % the markov_chain is exogenous
        markov_chains(loc).is_endogenous=false;
    end
    if isempty(blocks(paramBlock).listing(iparam).governing_chain)
        blocks(paramBlock).listing(iparam).governing_chain=const_loc;
    else
        loc=find(strcmp(blocks(paramBlock).listing(iparam).governing_chain,...
            markov_chain_names));
        blocks(paramBlock).listing(iparam).governing_chain=loc;
    end
end

% replace pseudofunctions in the model blocks
%--------------------------------------------
blocks(modelBlock).listing=parser.process_keywords(...
    blocks(modelBlock).listing,endovar_names);

%--------------------------------------------------------------------------

    function block=construct_list(block,rawline_,tokk)
        
        current_list={block.listing.name};
        
        rawline_without_description=parser.remove_description(rawline_);
        
        try
            
            end_game= contains(rawline_without_description,';');
            
        catch
            
            end_game=~isempty(strfind(rawline_without_description,';')); 
            
        end
        
        if end_game
            warning(['ending declaration listings with a'';'' is no longer required. In ',...
                file_name,' at line ',sprintf('%0.0f',line_number)])
        end
        if is_trigger
            if left_parenth_4_markov_chains_open
                error('new trigger cannot start before a listing a switching parameter block is complete')
            end
            [~,~,rawline_]=parser.look_around(tokk,rawline_);
            % reset the name of the markov chain
            current_markov_chain_name='';
            % check whether a parenthesis opens up
            test_left_par=strtok(rawline_);
            left_parenth_4_markov_chains_open= ~isempty(test_left_par) && ...
                strcmp(test_left_par(1),'(');
            if left_parenth_4_markov_chains_open
                new_markov_chain_tex_names={};
            end
        end
        while ~isempty(rawline_)
            if left_parenth_4_markov_chains_open && ~quote_active
                % close it if the next item is a ')'
                tokk=strtok(rawline_,strrep(DELIMITERS,')',''));
                if strcmp(tokk,')')
                    % make sure possibly only space was in between
                    right_par=find(rawline_==')',1,'first');
                    in_between=rawline_(1:right_par-1);
                    if isempty(in_between)||all(isspace(in_between))
                        % make sure it is not empty
                        if isempty(current_markov_chain_name)
                            error(['could not find the markov chain after "parameters" in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                        end
                        % make sure the number of tex names is adequate before
                        % pushing
                        %------------------------------------------------------
                        if ~isempty(new_markov_chain_tex_names)
                            if ~isequal(numel(new_markov_chain_tex_names),current_number_of_states)
                                error(['# descriptions for markov chain does not ',...
                                    'match the # of states ',file_name,' at line ',...
                                    sprintf('%0.0f',line_number)])
                            end
                            markov_chains(end).state_tex_names=new_markov_chain_tex_names;
                        end
                        % reset everything except the chain name as it will be
                        % used to dispatch the parameters
                        %------------------------------------------------------
                        current_number_of_states=[];
                        new_markov_chain_tex_names={};
                        left_parenth_4_markov_chains_open=false;
                        % update rawline_
                        %-----------------
                        rawline_=rawline_(right_par+1:end);
                    else
                        error(['"',in_between,'" not expected before closing of parenthesis in ',...
                            file_name,' at line ',sprintf('%0.0f',line_number)])
                    end
                end
            end
            [tokk,rest_]=strtok(rawline_,DELIMITERS);
            if ismember(tokk,blknames) &&  ~quote_active
                error([tokk,' not admissible as a declaration or block name ',...
                    'should start at the beginning of a line. In ',...
                    file_name,' at line ',sprintf('%0.0f',line_number)])
            end
            if ~isempty(tokk)
                if strcmp(tokk(1),'"')
                    if ~quote_active
                        quote_active=true;
                        tokk=tokk(2:end);
                    end
                end
                if quote_active
                    second_quote=strfind(tokk,'"');
                    if ~isempty(second_quote)
                        tokk=tokk(1:second_quote-1);
                    end
                    if isempty(tex_name)
                        tex_name=tokk;
                    else
                        if isempty(tokk)
                            if ~isempty(second_quote)
                                % we've just deleted a quote sign and so, look
                                % around a "
                                look_behind=parser.look_around('"',rawline_,true);
                            else
                                error([mfilename,':: parsing error. Please contact junior.maih@gmail.com'])
                            end
                        else
                            look_behind=parser.look_around(tokk,rawline_,true);
                        end
                        % preserve space for tex names
                        tex_name=[tex_name,look_behind,tokk]; %#ok<*AGROW>
                        tex_name=strtrim(tex_name);
                    end
                    if ~isempty(second_quote)
                        quote_active=false;
                        if left_parenth_4_markov_chains_open
                            new_markov_chain_tex_names=[new_markov_chain_tex_names,...
                                tex_name];
                        else
                            if ~isempty(block.listing(end).tex_name)
                                error([mfilename,':: found two consecutive tex names in ',...
                                    file_name,' at line ',sprintf('%0.0f',line_number)])
                            end
                            block.listing(end).tex_name=tex_name;
                        end
                        tex_name='';
                    end
                else
                    isv_name=isvarname(tokk);
                    if left_parenth_4_markov_chains_open
                        if isv_name
                            % markov chain name
                            %------------------
                            if isempty(current_markov_chain_name)
                                current_markov_chain_name=tokk;
                                % make sure it does not contain any _
                                if any(current_markov_chain_name=='_')
                                    error(['Markov chain names cannot contain "_" in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                                end
                                % make sure it is not the constant markov chain
                                if strcmp(current_markov_chain_name,'const')
                                    error(['"const" cannot be used as a name for a markov chain in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                                end
                            else
                                error([' two consecutive chain names in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                            end
                        else
                            % number of states
                            %------------------
                            tmp=str2double(tokk);
                            if isnan(tmp)||~(tmp==floor(tmp))
                                error([tokk,' is not a valid number of states in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                            end
                            if isempty(current_number_of_states)
                                current_number_of_states=tmp;
                                if isempty(current_markov_chain_name)
                                    error(['chain name must come before the number of states in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                                end
                                if current_number_of_states<=1
                                    error(['A declared markov chain must have more than one state. ',...
                                        'Remove the markov chain info if the parameters will ',...
                                        'remain unchanged across regimes. See ',file_name,...
                                        ' at line ',sprintf('%0.0f',line_number)])
                                end
                                % update the markov chain information
                                %-------------------------------------
                                loc_=find(strcmp(current_markov_chain_name,all_chain_names));
                                if isempty(loc_)
                                    markov_chains(end+1)=parser.initialize_markov_chain(current_markov_chain_name,current_number_of_states);
                                    all_chain_names={markov_chains.name};
                                else
                                    if ~isequal(markov_chains(loc_).number_of_states,current_number_of_states)
                                        error([current_markov_chain_name,' was previously declared to have ',...
                                            sprintf('%0.0f',markov_chains(loc_).number_of_states),' states but now has ',...
                                            sprintf('%0.0f',current_number_of_states),' in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                                    end
                                end
                            else
                                error([' two consecutive number of states in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                            end
                        end
                    else
                        if ~isv_name
                            error([' atom ''',tokk,''' in ',...
                                file_name,' at line ',sprintf('%0.0f',line_number),' is not a valid variable or parameter name'])
                        end
                        if exist([tokk,'.m'],'file')
                            disp([mfilename,':: (gentle warning): ',tokk,' is also a matlab function'])
                        end
                        [flag,mess]=parser.is_forbidden_name(tokk);
                        if flag
                            error([mess,' in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                        end
                        if ~isempty(current_list) &&  ismember(tokk,current_list)
                            error(['atom ''',tokk,''' has been declared twice in ',...
                                file_name,' at line ',sprintf('%0.0f',line_number)])
                        end
                        %----------- {name,tex_name,in-use-flag} --------------
                        nvars=numel(block.listing)+1;
                        block.listing(nvars)=parser.listing('name',tokk);
                        % update list in case many items are on the same line!!!
                        %-------------------------------------------------------
                        current_list={block.listing.name};
                        if strcmp(current_block_name,'parameters')
                            block.listing(nvars).governing_chain=current_markov_chain_name;
                            % exogenous transition probabilities can only be
                            % controlled by the const markov chain
                            %-----------------------------------------------
                            [istp,is_diagonal]=parser.is_transition_probability(tokk);
                            if is_diagonal
                                error(['Declaring diagonal transition probabilities is not allowed in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                            end
                            if isempty(current_markov_chain_name) % constant chain
                                block.listing(nvars).is_switching=false;
                            else
                                if istp
                                    error(['Transition probability "',tokk,'" cannot switch in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                                end
                                block.listing(nvars).is_switching=true;
                            end
                            block.listing(nvars).is_measurement_error=false;
                        end
                    end
                end
            end
            rawline_=rest_;
        end
    end
end

function [is_trigger,last_block_id]=check_block(tokk,trigger_map,last_block_id)
% {blknm,this_triggers{itrig},ii}
loc_=find(strcmp(tokk,trigger_map(:,2)));
is_trigger=false;
if ~isempty(loc_)
    is_trigger=true;
    last_block_id=trigger_map{loc_,3};
end
end

