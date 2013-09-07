function [transition_matrix_symbolic,transition_matrix_shadow,markov_chain_info,myifelseif]=...
    transition_probabilities(...
    input_list,parameters,markov_chains,shadow_tvp,shadow_definitions,probability_of_commitment)

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
parameter_names={parameters.name};
for i1=1:size(markov_chains,2)
    chain=markov_chains(i1).name;
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
                    if strcmp(chain,'loosecommit')
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
    transition_matrix_symbolic.(chain)=analytical_symbolic_form(NewQ,input_list,'symbolic');
    % list of symbols
    thisList=collect_symbolic_list(transition_matrix_symbolic.(chain),strcat(input_list,'_'));
    symb_list=union(symb_list,thisList);
end

[transition_matrix_symbolic,markov_chain_info,myifelseif]=parser.build_markov_regimes(transition_matrix_symbolic,markov_chains);

transition_matrix_shadow=analytical_symbolic_form(transition_matrix_symbolic,input_list,'analytic');
[nrows,ncols]=size(transition_matrix_shadow);
for irow=1:nrows
    for icol=1:ncols
        transition_matrix_shadow{irow,icol}=['Q(',sprintf('%0.0f',irow),',',sprintf('%0.0f',icol),')=',transition_matrix_shadow{irow,icol},';'];
    end
end
transition_matrix_shadow=[['Q=zeros(',sprintf('%0.0f',nrows),',',sprintf('%0.0f',ncols),');'];transition_matrix_shadow(:)];
transition_matrix_shadow=parser.code2func(transition_matrix_shadow,input_list(:)',{'Q'});
end
