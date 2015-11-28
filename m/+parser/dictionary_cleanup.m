function dictionary=dictionary_cleanup(dictionary,dynamic,static,...
    old_endo_names,logical_incidence)
cleanup_endogenous(old_endo_names,logical_incidence);

cleanup_exogenous();

cleanup_observables();

cleanup_model_types();

cleanup_parameters();

cleanup_equations(dynamic,static);

cleanup_definitions();

dictionary=orderfields(dictionary);

    function cleanup_endogenous(old_endo_names,logical_incidence)
        dictionary.endogenous=parser.greekify(dictionary.endogenous);
        endogenous=dictionary.endogenous;
        dictionary.endogenous=struct();
        dictionary.endogenous.name={endogenous.name};
        dictionary.endogenous.tex_name={endogenous.tex_name};
        dictionary.endogenous.current_name={endogenous.current_name};
        dictionary.endogenous.is_original=sparse(ismember(dictionary.endogenous.name,old_endo_names));
        dictionary.endogenous.number=numel(dictionary.endogenous.is_original);
        dictionary.endogenous.is_lagrange_multiplier=sparse(strncmp('MULT_',{endogenous.name},5));
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
        hbe={dictionary.parameters.name};
        hbe=hbe(strncmp(hbe,'hbe_param_',10));
        hbe=regexprep(hbe,'hbe_param_(\w+)','$1');
        dictionary.endogenous.is_hybrid_expect=ismember(dictionary.endogenous.name,hbe);
        
        % re-encode the auxiliary variables
        %------------------------------------
        dictionary.auxiliary_variables.sstate_solved=...
            union(dictionary.auxiliary_variables.model,dictionary.auxiliary_variables.ssmodel_solved);
        encode=@(x)vec(ismember(dictionary.endogenous.name,x)).';
        fields=fieldnames(dictionary.auxiliary_variables);
        for ifield=1:numel(fields)
            ff=fields{ifield};
            dictionary.auxiliary_variables.(ff)=encode(dictionary.auxiliary_variables.(ff));
        end
    end

    function cleanup_exogenous()
        dictionary.exogenous=parser.greekify(dictionary.exogenous);
        exogenous=dictionary.exogenous;
        dictionary.exogenous=struct();
        dictionary.exogenous.name={exogenous.name};
        dictionary.exogenous.tex_name={exogenous.tex_name};
        dictionary.exogenous.is_observed=sparse(ismember(dictionary.exogenous.name,{dictionary.observables.name}));
        dictionary.exogenous.number=full([sum(~dictionary.exogenous.is_observed),sum(dictionary.exogenous.is_observed)]);
        dictionary.exogenous.is_in_use=sparse([exogenous.is_in_use]);
        dictionary.exogenous.shock_horizon=sparse(dictionary.markov_chains.regimes_number,...
            sum(dictionary.exogenous.number));
    end

    function cleanup_observables()
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
    end

    function cleanup_model_types()
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
    end

    function cleanup_parameters()
        dictionary.parameters=parser.greekify(dictionary.parameters);
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
        dsge_var_id=strcmp(dictionary.parameters.name,parser.name4dsgevar());
        dictionary.is_dsge_var_model=any(dsge_var_id);
        dictionary.parameters.is_in_use(dsge_var_id)=true;
        dictionary.parameters.governing_chain=[parameters.governing_chain];
    end

    function cleanup_equations(dynamic,static)
        % dynamic
        dictionary.equations.dynamic=dynamic.model;
        dictionary.equations.shadow_dynamic=dynamic.shadow_model;
        dictionary.equations.number=numel(dynamic.model);
        % static
        dictionary.equations.static=static.model;
        dictionary.equations.shadow_static=static.shadow_model;
        dictionary.equations.shadow_steady_state_model=static.shadow_steady_state_model;
        dictionary.equations.shadow_steady_state_auxiliary_eqtns=static.shadow_steady_state_auxiliary_eqtns;
        dictionary.equations.shadow_fast_ssmodel=static.shadow_fast_ssmodel;
    end

    function cleanup_definitions()
        definitions=dictionary.definitions;
        dictionary.definitions=struct();
        dictionary.definitions.name={definitions.name};
        dictionary.definitions.dynamic={definitions.model};
        dictionary.definitions.dynamic=dictionary.definitions.dynamic(:);
        dictionary.definitions.shadow_dynamic={definitions.shadow};
        dictionary.definitions.shadow_dynamic=dictionary.definitions.shadow_dynamic(:);
        dictionary.definitions.number=numel(definitions);
    end
end
