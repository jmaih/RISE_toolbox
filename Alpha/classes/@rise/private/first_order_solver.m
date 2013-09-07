function [obj,retcode]=first_order_solver(obj,structural_matrices)
% solver options are 0 or 'msre_klein' (for constant-parameter models)
%                    1 or 'msre_gensys' (for constant-parameter models)
%                    2 or 'msre_aim' (for constant-parameter models)
%                    3 or 'functional_iteration'
%                    4 or 'newton_kronecker'
%                    5 or 'newton_system'
%                    6 or 'newton_kronecker_iteration'
%                    7 or 'loose_com' (for constant-parameter optimal policy models)
% the algorithm automatically switches to the constant parameter solver
% 0 (Klein) whenever possible

Defaults=struct();
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=Defaults;
    return
end
endo_nbr=obj.endogenous.number(2);
h=obj.markov_chains.regimes_number;
[Gplus,A0,Aminus,B,C]=msre_solvers.integrate_structure(structural_matrices,h); 
Q=obj.solution.Q; 
hidden_SS=obj.solution.ss;
hidden_BGP=obj.solution.bgp;
solve_expect_order=obj.options.solve_expect_order;

T0=[]; % endo_nbr x endo_nbr x h initial guess for T

% default solution for risk
gsig=repmat({zeros(endo_nbr,1)},1,h);

if obj.is_hybrid_expectations_model
    lambda=obj.parameter_values(obj.hybrid_expectations_lambda_id,:);
    [Gplus,Aminus]=hybrid_expectations(Gplus,Aminus,lambda);
end

if obj.is_sticky_information_model
    lambda=[obj.parameters(obj.sticky_information_lambda_id).startval];
    % the hidden steady state is computed in evaluate
    [Gplus,A0,Aminus,B,hidden_SS]=sticky_information(...
        Gplus,A0,Aminus,B,hidden_SS,lambda,...
        obj.forward_looking_ids);
    [Gplus,A0,Aminus,B,hidden_SS,hidden_BGP]=...
        sticky_information(Gplus,A0,Aminus,B,hidden_SS,hidden_BGP,lambda,obj.forward_looking_ids);
end

if obj.is_optimal_policy_model
    % do not take this steady state into consideration as
    % it automatically returns 0. the user might have given
    % a steady state file that returns different values.
    % Plus it does not account for balanced growth, which
    % is taken care of in the steady state computation.
    % Find the first non-zero diagonal probability
    dQ=diag(Q); first=find(dQ>0,1,'first');
    [TT,RR,~,retcode,obj.options]=dsge_lc_solve(Aminus{1},A0{1},...
        Gplus{first,first}/Q(first,first),B,...
        structural_matrices.planner.weights,...
        structural_matrices.planner.commitment,...
        structural_matrices.planner.discount,solve_expect_order,...
        obj.reordering_index,T0,obj.options);
    % change the name of the solver right here, right now!
    obj.options.solver='loose_com';
    if ~retcode
        [TT,RR,hidden_SS,hidden_BGP]=...
            loose_commitment_to_markov_switching(obj,TT,RR,hidden_SS,hidden_BGP);
    end
else
    sparam=obj.parameter_values(obj.parameters.is_switching,:);
    
    [TT,RR,gsig,theta_hat,retcode,obj.options] = msre_solve(Aminus,A0,Gplus,B,...
        C,sparam,obj.is_unique_steady_state,Q,T0,obj.options);
    
    if obj.is_sticky_information_model
        for st=1:h
            TT{st}=TT{st}(obj.reordering_index,obj.reordering_index);
            RR{st}=RR{st}(obj.reordering_index,:,:);
            hidden_SS{st}=hidden_SS{st}(obj.reordering_index,:);
            hidden_BGP{st}=hidden_BGP{st}(obj.reordering_index,:);
%             gsig{st}=gsig{st}(obj.reordering_index); % this is unnecessary since all are zero (no parameter switch), but for completeness
        end
    end
end

if ~retcode
    oldT=obj.solution.m_x;
    obj.solution.m_x=TT;
    if obj.options.check_stability && obj.markov_chains.regimes_number>1
        if (ischar(obj.options.solver)&& ...
                ~ismember(obj.options.solver,{'msre_gensys','msre_aim'}))||...
                (isnumeric(obj.options.solver)&& ...
                ~ismember(obj.options.solver,[1,2]))
            if ~obj.is_stable_system
                obj.solution.m_x=oldT;
                retcode=25; % system unstable
            end
        end
    end
    if ~retcode
        % This is where I should apply the shock restrictions before passing everything
        % to obj.solution.m_e... The format of R is endo_nbr x exo_nbr x solve_expect_order x regime
        if solve_expect_order>1 && ~isempty(obj.options.shock_properties)
            nshocks=sum(obj.exogenous.number);
            for ii=1:nshocks
                shock=obj.options.shock_properties(ii).name;
                shock_id=strcmp(shock,{obj.varexo.name});
                horizon=obj.options.shock_properties(ii).horizon;
                for istate=1:numel(RR)
                    RR{istate}(:,shock_id,horizon+1:end)=0;
                end
            end
        end
        obj.solution.m_e=RR;
        obj.solution.ss=hidden_SS;
        obj.solution.bgp=hidden_BGP;
        obj.solution.m_sig=gsig;
        if ~obj.is_optimal_policy_model
            % this will be used in higher-order approximations
            obj.solution.theta_hat=theta_hat;
        end
        % store the following information so that we don't have to
        % resolve the model if we don't have to. I probably should put
        % an extra option called current_solver so that the model is
        % resolved if the solver changes.
        obj.current_solution_parameters=obj.parameter_values;
    end
end

    function [Aplus_he,Aminus_he]=hybrid_expectations(Gplus,Aminus,lambda)
        Aminus_he=Aminus;
        Aplus_he=Gplus;
        for s0=1:h
            Aminus_he{s0}=Aminus{s0};
            for s1=1:h
                Aminus_he{s0}=Aminus_he{s0}+(1-lambda(s1))*Gplus{s0,s1};
                Aplus_he{s0,s1}=lambda(s1)*Gplus{s0,s1};
            end
        end
    end

    function [Aplus_si,A0_si,Aminus_si,B_si,SS_si,BGP_si]=...
            sticky_information(Gplus,A0,Aminus,B,SS,BGP,lambda,truly_forward_id)
        [old_endo_nbr,exo_nbr]=size(B{1});
        new_endo_nbr=endo_nbr+numel(truly_forward_id);
        % initialize output matrices
        Aplus_si=Gplus;
        A0_si=A0;
        Aminus_si=Aminus;
        B_si=B;
        SS_si=SS;
        BGP_si=BGP;
        % define the steady_state and balance growth paths locations
        fw=numel(truly_forward_id);
        
        % the variables are arranged as : [Y',SI']' with the sticky information
        % variables below, Y_{t+1}=lamb*Y_{t+1}+(1-lamb)*SI_{t-1} and
        % SI_t=E_t(Y_{t+1})
        
        for s0=1:h
            A0_si{s0}=eye(new_endo_nbr);
            A0_si{s0}(1:old_endo_nbr,1:old_endo_nbr)=A0{s0};
            Aminus_si{s0}=zeros(new_endo_nbr);
            Aminus_si{s0}(1:old_endo_nbr,1:old_endo_nbr)=Aminus{s0};
            B_si{s0}=[B{s0};zeros(fw,exo_nbr)];
            for s1=1:h
                Aminus_si{s0}(1:old_endo_nbr,old_endo_nbr+1:end)=...
                    Aminus_si{s0}(1:old_endo_nbr,old_endo_nbr+1:end)+(1-lambda(s1))*Gplus{s0,s1}(:,truly_forward_id);
                Aplus_si{s0,s1}(1:old_endo_nbr,1:old_endo_nbr)=lambda(s1)*Gplus{s0,s1};
                Aplus_si{s0,s1}(1:old_endo_nbr,truly_forward_id)=-eye(fw);
            end
            SS_si{s0}=[SS{s0};SS{s0}(truly_forward_id)];
            BGP_si{s0}=[BGP{s0};BGP{s0}(truly_forward_id)];
        end
    end
end

