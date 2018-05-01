function [block,dictionary,static]=capture_equations(dictionary,listing,block_name,static)
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

if nargin<4
    static=struct();
end
function_on=false;
time_on=false;
last_status='';
fill_time='';
time_opening='';
def_flag=false;
endo_switch_flag=false;
mcp_flag=false;
equation=initialize_equation();%cell(2,0);
block=cell(0,5);
DELIMITERS=parser.delimiters();

chain_names={dictionary.markov_chains.name};
chain_states_number=[dictionary.markov_chains.number_of_states];
parameter_names={dictionary.parameters.name};
parameter_govern=[dictionary.parameters.governing_chain];

% profile off
% profile on
% [listing,nlist]=parser.insert_definitions(listing, dictionary.definitions_inserted);
% profile off
% profile viewer
% keyboard
nlist=size(listing,1);

swp=switching_status_of_parameters(dictionary);

nblks=0;
blocks_to_discard_coz_they_are_defs=false(nblks,1);
definitions_loc=struct();

for ii=1:nlist
    [block,equation]=capture_equations_engine(block,listing(ii,:),block_name,equation);
end

if dictionary.definitions_inserted
    block=block(~blocks_to_discard_coz_they_are_defs,:);
    % remove all definitions from the list
    dictionary.definitions={};
end

    function eqtn=initialize_equation()
        eqtn=struct('max_lag',0,'max_lead',0,...
            'eqtn',{cell(2,0)},...
            'type','normal',...
            'is_def',false,... % definitions
            'is_tvp',false,... % endogenous probabilities
            'is_mcp',false); % complementarity condition
    end
    function [block,equation]=capture_equations_engine(block,cell_info,block_name,equation)
        % if a parameter has a lead, create an auxiliary variable for it
        % if a variable has a lead or a lag greater than one, replace it
        % with the corresponding auxiliary variable and update the lead or
        % lag of the variable in the list and then create auxiliary
        % equations when the parsing of the model is done.
        
        % store this useful information separately
        %-----------------------------------------
        reset_lists();
        
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
            
            if isempty(equation.eqtn)
                endo_switch_flag=false;
                def_flag=false;
                if time_on
                    error([mfilename,':: new equation starting without finished time index in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                end
                if function_on
                    error([mfilename,':: new equation starting with earlier function not closed in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
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
                                error([mfilename,':: attributes to ',...
                                    block_name,' block no longer permitted in file ',...
                                    file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                        end
                        break % exit while ~isempty(rest_) loop
                    end
                end
                
                if ~isempty(tokk)
                    tok_status=dictionary.determine_status(tokk,dictionary);
                    if strcmp(tok_status,'param')
                        if strcmp(block_name,'exogenous_definition')
                            error([mfilename,':: exogenous definitions cannot contain parameters in file ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                        position=strcmp(tokk,{dictionary.parameters.name});
                        dictionary.parameters(position).is_in_use=true;
                    elseif strcmp(tok_status,'x')
                        position=strcmp(tokk,{dictionary.exogenous.name});
                        dictionary.exogenous(position).is_in_use=true;
                    end
                    is_lhs_def=false;
                    if (strcmp(tok_status,'#') && isempty(equation.eqtn))||...
                            strcmp(tok_status,'!')||...
                            strcmp(tok_status,'?')
                        def_flag=strcmp(tok_status,'#');
                        if strcmp(block_name,'exogenous_definition')
                            error([mfilename,':: the exogenous definition block cannot contain ''#'' or ''!'' or ''?'' ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                        endo_switch_flag=strcmp(tok_status,'!');
                        mcp_flag=strcmp(tok_status,'?');
                        rest_=rest1;
                        [tokk,rest1]=strtok(rest_,DELIMITERS); %#ok<*STTOK>
                        tok_status=dictionary.determine_status(tokk,dictionary);
                        if ~mcp_flag && ~strcmp(tok_status,'unknown')
                            if strcmp(tok_status,'f')
                                disp([mfilename,':: (gentle warning): ',tokk,' is also a matlab function'])
                            else
                                error([mfilename,':: string ''',tokk,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_),' cannot have multiple types'])
                            end
                        end
                        if def_flag
                            dictionary.definitions=[dictionary.definitions;{tokk}];
                            is_lhs_def=true;
                            equation.is_def=true;
                            equation.type='def';
                            if isempty(tokk)
                                error([mfilename,':: definition missing in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            definitions_loc.(tokk)=nblks+1;
                        elseif mcp_flag
                            equation.is_mcp=true;
                            equation.type='mcp';
                            % everything shall be parse as normal
                            % equations. but still, perhaps this is the
                            % place where to stamp the equations for their
                            % type=
                        elseif endo_switch_flag
                            equation.is_tvp=true;
                            equation.type='tvp';
                            [istp,isdiagonal,chain_name]=parser.is_transition_probability(tokk);
                            if ~istp
                                error([mfilename,':: string ''',tokk,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_),' is not an appropriate name for an endogenous switching probability'])
                            end
                            if isdiagonal
                                error([mfilename,':: "',tokk,'" is a diagonal transition probabilitiy. Only off-diagonal elements are allowed. Check ',file_name_,' at line ',sprintf('%0.0f',iline_),' is not an appropriate name for an endogenous switching probability'])
                            end
                            dictionary.time_varying_probabilities=[dictionary.time_varying_probabilities,{tokk}];
                            ch_names={dictionary.markov_chains.name};
                            loc_chain=find(strcmp(chain_name,ch_names));
                            if isempty(loc_chain)
                                error([mfilename,':: markov chain "',chain_name,'" has not be declared. In ',file_name_,' at line ',sprintf('%0.0f',iline_),' is not an appropriate name for an endogenous switching probability'])
                            end
                            chain_status=dictionary.markov_chains(loc_chain).is_endogenous;
                            if isnan(chain_status)
                                dictionary.markov_chains(loc_chain).is_endogenous=true;
                            elseif ~isequal(chain_status,true)
                                error([mfilename,':: markov chain "',chain_name,'" was previously found to be exogenous and now is endogenous. In ',file_name_,' at line ',sprintf('%0.0f',iline_),' is not an appropriate name for an endogenous switching probability'])
                            end
                        else
                            error([mfilename,':: parsing error in ',file_name_,' at line ',sprintf('%0.0f',iline_),' please report this to junior.maih@gmail.com'])
                        end
                        % update the status above
                        tok_status=dictionary.determine_status(tokk,dictionary);
                        if ~isempty(equation.eqtn)
                            error([mfilename,':: # and ! can only occur at the beginning of an equation check in ',file_name_,' line ',sprintf('%0.0f',iline_)])
                        end
                    elseif strcmp(tok_status,'unknown')
                        if parser.is_transition_probability(tokk) && strcmp(block_name,'model')
                            error([mfilename,':: equations for endogenous transition probabilities must start with a "!" in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        else
                            error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                    end
                    
                    if def_flag && (strcmp(tok_status,'y')||strcmp(tok_status,'x')) && ~dictionary.definitions_inserted
                        error([mfilename,':: definitions cannot contain variables. ',...
                            'If this is an endogenous parameter, declare the definitions as ',...
                            'parameters and use a steady state file. check ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                    end
                    
                    if ~endo_switch_flag && strcmp(tok_status,'tvp')
                        error([mfilename,':: model equations cannot contain endogenous switching probabilities. check ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                    end
                    left_operator=parser.look_around(tokk,rest_);
                    
                    for i1=2:length(left_operator)
                        first=dictionary.determine_status(left_operator(i1-1),dictionary);
                        if strcmp(first,'unknown')
                            error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                        second=dictionary.determine_status(left_operator(i1),dictionary);
                        if strcmp(second,'unknown')
                            error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                        check_validity([first,second],file_name_,iline_,block_name);
                    end
                    if ~isempty(left_operator)
                        if ~isempty(tokk)
                            third=dictionary.determine_status(left_operator(end),dictionary);
                            if strcmp(third,'unknown')
                                error([mfilename,':: unknown string ''',tokk,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            check_validity([third,tok_status],file_name_,iline_,block_name)
                        end
                        % always deal with time before dealing with
                        % functions as time potentially removes some
                        % parentheses and the functions do not.
                        if (...
                                (strcmp(last_status,'y')||...
                                strcmp(last_status,'x')||...
                                (strcmp(last_status,'param')&& ~strcmp(block_name,'parameter_restrictions'))...
                                ) &&...
                                (strcmp(left_operator(1),'(')||strcmp(left_operator(1),'{')))% check the time
                            time_on=true;
                            time_opening=left_operator(1);
                            fill_time=[fill_time,left_operator(2:end)];
                            left_operator='';
                        elseif time_on && (strcmp(left_operator(1),')')||strcmp(left_operator(1),'}'))
                            if isempty(fill_time)
                                error([mfilename,':: missing lead or lag in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            equation=update_leads_lags(equation,fill_time);
                            fill_time='';
                            time_on=false;
                            time_closing=left_operator(1);
                            if ~((strcmp(time_opening,'(') && strcmp(time_closing,')'))||...
                                    (strcmp(time_opening,'{') && strcmp(time_closing,'}')))
                                error([mfilename,':: time opened with ''',time_opening,''' and closed with ''',time_closing,''' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
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
                            error([mfilename,':: parenthesis mismatch in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                        
                        if ~isempty(left_operator)
                            equation.eqtn=[equation.eqtn,{left_operator,[]}']; %#ok<*AGROW>
                        end
                    end
                    if time_on && strcmp(tok_status,'n')
                        fill_time=[fill_time,tokk];
                    else
                        % if it is a definition, get rid of it if the
                        % definitions are to be inserted
                        %----------------------------------------------
                        if ~is_lhs_def && strcmp(tok_status,'def') && dictionary.definitions_inserted
                            old_eqtn=block{definitions_loc.(tokk),1};
                            middle_man=old_eqtn;
                            % remove the first part
                            %----------------------
                            middle_man(:,1)=[];
                            % remove the equality sign
                            %-------------------------
                            middle_man{1,1}(1)=[];
                            if isempty(middle_man{1,1})
                                middle_man(:,1)=[];
                            end
                            % remove the semicolon
                            %---------------------
                            middle_man{1,end}(end)=[];
                            if isempty(middle_man{1,end})
                                middle_man(:,end)=[];
                            end
                            block_def=[{'(',[]}',middle_man,{')',[]}'];
                            equation.eqtn=[equation.eqtn,block_def];
                        else
                            equation.eqtn=[equation.eqtn,{tokk,[]}'];
                        end
                        if (strcmp(tok_status,'y')||strcmp(tok_status,'x')||...
                                (strcmp(tok_status,'param') && swp.(tokk))...
                                )
                            equation.eqtn{2,end}=0;
                        end
                        fill_time='';
                    end
                    last_status=tok_status;
                    rest_=rest1;
                else
                    if time_on && (strcmp(rest_(1),')')||strcmp(rest_(1),'}'))
                        rest_=rest_(2:end);
                        %====================
                        equation=update_leads_lags(equation,fill_time);
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
                        rest_(isspace(rest_))=[];
                        if ~isempty(rest_)
                            equation.eqtn=[equation.eqtn,{rest_,[]}'];
                        end
                        rest_='';
                    end
                    last_status=dictionary.determine_status(equation.eqtn{1,end}(end),dictionary);
                end
            end
            if ~isempty(equation.eqtn)
                if strcmp(equation.eqtn{1,end}(end),';')
                    % we've reach the end of the equation, validate it,
                    % load it and reinitialize.
                    equation.eqtn=validate_equation(equation.eqtn,...
                        max(abs([equation.max_lag,equation.max_lead])));
                    fast_ss_form=[];
                    pound_key=find(strcmp(equation.eqtn(1,:),'#'));
                    if ~isempty(pound_key)
                        fast_ss_form=equation.eqtn(:,pound_key+1:end);
                        equation.eqtn=equation.eqtn(:,1:pound_key-1);
                    end
                    block=[block;...
                        {equation.eqtn,equation.max_lag,equation.max_lead,equation.type,fast_ss_form}];
                    nblks=nblks+1;
                    if equation.is_def && dictionary.definitions_inserted
                        blocks_to_discard_coz_they_are_defs(nblks,1)=true;
                    else
                        blocks_to_discard_coz_they_are_defs(nblks,1)=false;
                    end
                    equation=initialize_equation();
                    % clear the mcp_flag
                    %--------------------
                    if mcp_flag
                        mcp_flag=false;
                    end
                end
            end
            % make sure the list is up to date for status determination
            %-----------------------------------------------------------
            reset_lists()
        end
        function reset_lists()
            do_one_list('endogenous')
            do_one_list('parameters')
            do_one_list('exogenous')
            function do_one_list(type_)
                fake_name=[type_,'_list'];
                if ~isfield(dictionary,fake_name)||...
                        numel(dictionary.(fake_name))~=numel(dictionary.(type_))
                    dictionary.(fake_name)={dictionary.(type_).name};
                end
            end
        end
        function equation=validate_equation(equation,max_lead_lag)
            % if there is a # in the equation, split the equation in to two
            % validate each part and reconnect
            pound_key_=find(strcmp(equation(1,:),'#'));
            if ~isempty(pound_key_)
                if ~strcmp(block_name,'model')
                    error(['# sign only allowed in model block in ',...
                        file_name_,' at line ',sprintf('%0.0f',iline_)])
                elseif numel(pound_key_)>1
                    error([' Only one # sign allowed per equation. In ',...
                        file_name_,' at line ',sprintf('%0.0f',iline_)])
                elseif mcp_flag||endo_switch_flag||def_flag
                    error([' Complementarities, endogenous switching probs ',...
                        ' and definitions cannot contain a # sign. In ',...
                        file_name_,' at line ',sprintf('%0.0f',iline_)])
                end
                % add a semicolon to the first part and validate
                first_part=[equation(:,1:pound_key_-1),{';';[]}];
                second_part=equation(:,pound_key_+1:end);
                % validate the bgp/sstate part
                first_part=validate_equation(first_part,max_lead_lag);
                second_part=validate_equation(second_part,max_lead_lag);
                % reconnect both
                equation=[first_part,{'#';[]},second_part];
            end
            if mcp_flag
                if max_lead_lag
                    % 3- cannot contain lags or leads
                    error(['A complementarity constraint cannot contain leads or lags in ',...
                        file_name_,' at line ',sprintf('%0.0f',iline_)])
                end
                mcp_test_passed=false;
                mcp_inequalities={'>=','<='...,'>','<','ge','gt','le','lt',
                    };
            end
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
                if mcp_flag
                    % Restrictions for complementarity constraints
                    % 1- no switching parameter
                    loc=strcmp(equation{1,ic},{dictionary.parameters.name});
                    if any(loc) && dictionary.parameters(loc).is_switching
                        error(['A complementarity constraint cannot contain switching parameters in ',...
                            file_name_,' at line ',sprintf('%0.0f',iline_)])
                    end
                    % 2- one of the following ge gt le lt
                    for iii=1:numel(mcp_inequalities)
                        loc=strfind(equation{1,ic},mcp_inequalities{iii});
                        if ~isempty(loc)
                            break
                        end
                    end
                    if ~isempty(loc)
                        if ~mcp_test_passed
                            if (ic==1||ic==size(equation,2)-1)
                                error(['inequality mis-placed in complementarity constraint in ',...
                                    file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            mcp_test_passed=true;
                        else
                            error(['A complementarity constraint cannot contain more than one inequality in ',...
                                file_name_,' at line ',sprintf('%0.0f',iline_)])
                        end
                    end
                elseif strcmp(block_name,'parameter_restrictions')
                    % Restrictions for the parameter restriction block
                    % 1- no variable
                    if ismember(equation{1,ic},{dictionary.endogenous.name})
                        error([mfilename,':: no variable allowed in the parameter_restrictions block'])
                    end
                    % 2- if parenthesis after a parameter, then one of the elements
                    % inside the parenthesis must be a chain name and the other a
                    % numeric.
                    locs=find(strcmp(equation{1,ic},parameter_names));
                    if ~isempty(locs) % chain_states_number={dictionary.markov_chains.number_of_states};
                        nom_de_chaine=chain_names{parameter_govern(locs)};
                        etats_de_la_chaine=1:chain_states_number(parameter_govern(locs));
                        state_=1;
                        chain_='const';
                        if  ic<size(equation,2) && strcmp(equation{1,ic+1}(1),'(')
                            if strcmp(nom_de_chaine,'const')
                                error(['parameter ',equation{1,ic},' is controlled by markov chain "const". So an opening parenthesis is unexpected in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            chain_=equation{1,ic+2};
                            if ~strcmp(chain_,nom_de_chaine) % chain
                                error(['parameter ',equation{1,ic},' is not controlled by a markov chain with name ',chain_,' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            if ~strcmp(equation{1,ic+3},',') % remove comma
                                error([mfilename,':: comma missing between chain and state for parameter ',equation{1,ic},' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            state_=eval(equation{1,ic+4});
                            if isnan(state_)||~(isnumeric(state_) && ismember(state_,etats_de_la_chaine))
                                error([mfilename,':: state invalid for parameter ',equation{1,ic},' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            if ~strcmp(equation{1,ic+5}(1),')') % closing parenthesis
                                error([mfilename,':: closing parenthesis missing for parameter ',equation{1,ic},' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                            end
                            % remove all processed elements
                            equation{1,ic+1}(1)=[]; % left parenthesis
                            equation{1,ic+2}=''; % state or chain
                            equation{1,ic+3}=''; % comma
                            equation{1,ic+4}=''; % state or chain
                            equation{1,ic+5}(1)=[]; % right parenthesis
                        end
                        % we carefully record the chain and the state for
                        % later processing
                        equation{2,ic}={chain_,state_};
                    end
                    % 3- parameter must be effectively controlled by the
                    % chain in the parameterization block
                end
            end
            % an equation cannot have multiple equality signs
            if size(equality_signs,2)>1
                error([mfilename,':: multiple equality signs found in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
            elseif size(equality_signs,2)==1
                    % validate the lhs and the rhs of the equality separately
                    lhs=equation(:,1:equality_signs{1}-1);
                    rhs=equation(:,equality_signs{1}+1:end);
                    if isempty(rhs)||isempty(lhs)
                        error([mfilename,':: left-hand side or right-hand side badly specified in ',...
                            file_name_,' at line ',sprintf('%0.0f',iline_)])
                    end
                    % split the middle cell and isolate the equality sign
                    string=equation{1,equality_signs{1}};
                    left_string=[];
                    right_string=[];
                    if mcp_flag
                        % locate > or < and =
                        inequal=strfind(string,'<');
                        inequal_type=1;
                        if isempty(inequal)
                            inequal=strfind(string,'>');
                            inequal_type=2;
                        end
                        equal_=strfind(string,'=');
                        lhs=[lhs,transpose({string(1:inequal-1),[]})];
                        rhs=[transpose({string(equal_+1:end),[]}),rhs];
                        if inequal_type==1
                            lhs=[transpose({'-(',[]}),lhs,transpose({')',[]})];
                        end
                    else
                        left_string=string(1:equality_signs{2}-1);
                        right_string=string(equality_signs{2}+1:end);
                    end
                    if ~isempty(left_string)
                        % add to the lhs
                        lhs=[lhs,transpose({left_string,[]})];
                    end
                    msg=count_parentheses(lhs(1,:));
                    if ~isempty(msg)
                        error([mfilename,':: ',msg,' on left hand side in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                    end
                    
                    if ~isempty(right_string)
                        % add to the rhs but before !
                        rhs=[transpose({right_string,[]}),rhs];
                    end
                    msg=count_parentheses(rhs(1,:));
                    if ~isempty(msg)
                        error([mfilename,':: ',msg,' on right hand side in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
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
                    if mcp_flag
                        equation{1,end}=[equation{1,end}(1:end-1),'>=0;'];
                    end
            end
            if ~checked
                % then there was no equality sign and so, check the whole
                % equation in one go
                msg=count_parentheses(equation(1,:));
                if ~isempty(msg)
                    error([mfilename,':: ',msg,' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
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
                    error([mfilename,':: ',vnamex,' must be exogenous and observed in ',file_name_,' at line ',sprintf('%0.0f',iline_)'])
                end
                % It must be dated at time 0
                if ~isequal(equation{2,1},0)
                    error([mfilename,':: ',vnamex,' must be dated at time 0 in ',file_name_,' at line ',sprintf('%0.0f',iline_)'])
                end
                % the second element must be an equality sign
                if ~strcmp(equation{1,2}(1),'=')
                    error([mfilename,':: the second token must be a ''=''  in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
                end
                % no rhs variable can be dated at a time beyond 0
                for icol=3:size(equation,2)
                    if ~isempty(equation{2,icol}) && equation{2,icol}>0
                        error([mfilename,':: right-hand side variables cannot be dated in the future in ',file_name_,' at line ',sprintf('%0.0f',iline_)'])
                    end
                end
            end
            if mcp_flag
                if ~mcp_test_passed
                    disp(mcp_inequalities)
                    error([' A complementarity condition must contain one of the above ''=''  in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
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
        function equation=update_leads_lags(equation,leadorlag)
            leadorlag=eval(leadorlag);
            if ~isnumeric(leadorlag)||~isequal(leadorlag,floor(leadorlag))
                error([mfilename,':: time syntax error in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
            end
            equation.eqtn{2,end}=leadorlag;
            equation.max_lag=min(equation.max_lag,leadorlag);
            equation.max_lead=max(equation.max_lead,leadorlag);
        end
    end
    function check_validity(syntax,file_name_,iline_,block_name)
        special_syntax=any(strcmp(syntax,dictionary.syntax_special));
        good=any(strcmp(syntax,dictionary.syntax_typical))||...
            (any(strcmp(syntax,dictionary.syntax_time)) && ~strcmp(block_name,'steady_state_model'))||...
            (any(strcmp(syntax,dictionary.syntax_function)) && function_on)||...
            (special_syntax && strcmp(block_name,'steady_state_model'))||...
            (special_syntax && strcmp(block_name,'parameter_restrictions'))||...
            (special_syntax && strcmp(block_name,'model') && mcp_flag);
        if ~good
            error([mfilename,':: wrong syntax ',syntax,' in ',file_name_,' at line ',sprintf('%0.0f',iline_)])
        end
    end
end

function p=switching_status_of_parameters(dic)
% control parameters
%-------------------
p=struct();
for ipar=1:numel(dic.parameters)
    
    p.(dic.parameters(ipar).name)=dic.parameters(ipar).is_switching;
    
end
end
