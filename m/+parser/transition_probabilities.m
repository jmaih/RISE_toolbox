function [transition_matrix_symbolic,transition_matrix_shadow,...
    markov_chain_info,myifelseif,sig_pos]=transition_probabilities(...
input_list,parameter_names,markov_chains,shadow_tvp,...
shadow_definitions,probability_of_commitment)
% INTERNAL FUNCTION
%

% sig_pos : location of the perturbation parameter in the parameter list:
% we differentiate with respect to the perturbation parameter since it can
% appear directly in the structural model

if ~isempty(shadow_definitions)
    
    shadow_definitions={shadow_definitions};
    
end

shadow_tvp_left_right=regexp(shadow_tvp,'=','split');

shadow_tvp_left_right=vertcat(shadow_tvp_left_right{:});

% remove semicolon
if ~isempty(shadow_tvp_left_right)
    
    shadow_tvp_left_right=strrep(shadow_tvp_left_right,';','');

end

transition_matrix_symbolic=struct();

symb_list={};

chain_names={markov_chains.name};

nchains=numel(chain_names);

for i1=1:nchains

    chain=chain_names{i1};
    
    states_nbr_ii=markov_chains(i1).number_of_states;
    
    endogeneity=markov_chains(i1).is_endogenous;
    
    if isnan(endogeneity)
    
        error(['endogeneity/exogeneity status of markov chain ',chain,' could not be determined. Most likely, the transition probs remain undeclared'])
    
    end
    
    exo_flag=~endogeneity;
    
    NewQ=cell(states_nbr_ii);
    
    for s1=1:states_nbr_ii
    
        cumul='0';
        
        for s2=1:states_nbr_ii
        
            if s1~=s2
            
                prob_name=[chain,'_tp_',sprintf('%0.0f',s1),'_',sprintf('%0.0f',s2)];
                
                if exo_flag
                
                    if strcmp(chain,parser.loose_commit())
                    
                        if s1==1 && s2==2
                        
                            new_tp=['(1-(',probability_of_commitment,'))'];
                        
                        else
                            
                            new_tp=['(',probability_of_commitment,')'];
                        
                        end
                        
                    else
                        
                        loc=find(strcmp(prob_name,parameter_names));
                        
                        if isempty(loc)
                        
                            error([mfilename,':: exogenous transition probability ',...
                                prob_name,' uncharacterized'])
                        
                        end
                        
                        new_tp=['param(',sprintf('%0.0f',loc),')'];
                    
                    end
                    
                else
                    
                    prob_loc= find(strcmp(shadow_tvp_left_right(:,1),prob_name));
                    
                    if isempty(prob_loc)
                    
                        error([mfilename,':: endogenous transition probability ',...
                            prob_name,' uncharacterized'])
                    
                    end
                    
                    new_tp=shadow_tvp_left_right{prob_loc,2};
                
                end
                
                NewQ{s1,s2}=new_tp;
                
                cumul=strcat(cumul,'+',new_tp);
            
            end
            
        end
        
        if strcmp(cumul,'0')
        
            NewQ{s1,s1}='1';
        
        else
            
            NewQ{s1,s1}=['1-(',cumul(3:end),')'];
        
        end
        
    end
    % those probabilities could be functions of
    % definitions. so take care of that right here.
    NewQ=parser.replace_definitions(NewQ,shadow_definitions);
    
    transition_matrix_symbolic.(chain)=parser.analytical_symbolic_form(NewQ,input_list,'symbolic');
    % list of symbols
    thisList=parser.collect_symbolic_list(transition_matrix_symbolic.(chain),strcat(input_list,'_'));
    
    symb_list=union(symb_list,thisList);
    
end

others=parser.perturbation_control_param_names();

pos=locate_variables(others,parameter_names,true);

sip=sum(isnan(pos));

% do not allow nan to pass otherwise the numerical differentiator will
% return wrong results
if sip==numel(pos)
    
    coef=[];
    
    sig_pos=[];
    
elseif sip==0
    
    panal=parser.parameters_analytical_form(pos);
    
    sig_pos=pos(1);
    
    sig=panal{1};
    
    iota1=panal{2};
    
    iota2=panal{3};
    
    coef=['(1+',iota1,'*',iota2,'*(',sig,'-1))'];
    
    coef=parser.symbolic_model(coef,'param');
    
else
    
    error('Unpleasant arithmetics: this cannot happen')

end

[~,markov_chain_info,myifelseif]=parser.build_markov_regimes(transition_matrix_symbolic,...
markov_chains,true,coef);

transition_matrix_shadow=struct();

for i1=1:nchains
    
    chain=chain_names{i1};
    
    cloc=strcmp(chain,markov_chain_info.regimes(1,:));
    
    regvec=cell2mat(markov_chain_info.regimes(2:end,cloc));
    
    transition_matrix_shadow.(chain)=transtion_matrix_to_function(...
        transition_matrix_symbolic.(chain),regvec);
    
end

    function Qfunc=transtion_matrix_to_function(TM,regimevec)

        Qfunc=['@(',cell2mat(strcat(input_list(:)',','))];
        
        Qfunc=[Qfunc(1:end-1),')['];
        
        for irow_=1:size(TM,1)
            
            state=irow_;
            
            ping=num2cell(find(regimevec==state));
            
            ping=cellfun(@int2str,ping(:)','uniformOutput',false);
            
            if numel(ping)>1
                
                ping=cell2mat(strcat(ping,','));
                
                ping=['[',ping(1:end-1),']'];
                
            else
                
                ping=ping{1};
                
            end
        
            tmp=cell2mat(strcat(TM(state,:),','));
            
            tmp=regexprep(tmp,'(def|y|param|ss|x)_(\d+)',['lp($1($2,:),',ping,')']);
            
            Qfunc=[Qfunc,tmp(1:end-1),';']; %#ok<AGROW>
        
        end
        
        Qfunc=[Qfunc(1:end-1),']'];
        
        Qfunc=str2func(Qfunc);
    
    end

end
