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
endo_nbr=obj.endogenous.number(end);
h=obj.markov_chains.regimes_number;
[Gplus,A0,Aminus,B,C]=msre_solvers.integrate_structure(structural_matrices,h); 
Q=obj.solution.Q; 
hidden_SS=obj.solution.ss;
hidden_BGP=obj.solution.bgp;
shock_horizon=max(obj.exogenous.shock_horizon);

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
        structural_matrices.planner.discount,shock_horizon,...
        obj.reordering_index,T0,obj.options);
    % change the name of the solver right here, right now!
    obj.options.solver='loose_com';
    if ~retcode
        [TT,RR,hidden_SS,hidden_BGP]=...
            loose_commitment_to_markov_switching(obj,TT,RR,hidden_SS,hidden_BGP);
    end
else
    sparam=obj.parameter_values(obj.parameters.is_switching,:);
	shock_horizon=max(obj.exogenous.shock_horizon);
    if obj.options.solve_accelerate
        [TT,RR,gsig,theta_hat,retcode,obj.options] = solve_small_msre_system(Aminus,A0,Gplus,B,...
            C,sparam,obj.is_unique_steady_state,Q,T0,shock_horizon,obj.options,...
            obj.endogenous.is_static);
    else
        [TT,RR,gsig,theta_hat,retcode,obj.options] = msre_solve(Aminus,A0,Gplus,B,...
            C,sparam,obj.is_unique_steady_state,Q,T0,shock_horizon,obj.options);
    end
    
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
        % to obj.solution.m_e... The format of R is endo_nbr x exo_nbr x shock_horizon x regime
        if shock_horizon>1
            nshocks=sum(obj.exogenous.number);
            for ii=1:nshocks
                horizon=obj.exogenous.shock_horizon(ii);
                for istate=1:numel(RR)
                    RR{istate}(:,ii,horizon+1:end)=0;
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

function [m_x,m_e,m_sig,theta_hat,retcode,options]=solve_small_msre_system(...
    Aminus,A0,Gplus,B,C,sparam,is_unique_steady_state,Q,T0,shock_horizon,options,stat_cols)
% solve a smaller system, partialling out static variables and then
% derive the solution for the bigger
nstates=numel(A0);
[AA0,AAminus,GBARplus,BB,CC]=msre_computational_savings();
endo_nbr=numel(stat_cols);
orig_order=1:endo_nbr;
dyn_cols=~stat_cols;
static_endo_nbr=sum(stat_cols);
dd=static_endo_nbr+1:endo_nbr;
if ~isempty(T0)
    T0=T0(dyn_cols,dyn_cols,:);
end
AA0_=AA0;AAminus_=AAminus;GBARplus_=GBARplus;BB_=BB;CC_=CC;
for istate=1:nstates
    AA0_{istate}=AA0{istate}(dd,dd);
    AAminus_{istate}=AAminus_{istate}(dd,dd);
    BB_{istate}=BB_{istate}(dd,:,:);
    for jstate=1:nstates
        CC_{istate,jstate}=CC_{istate,jstate}(dd,:);
        GBARplus_{istate,jstate}=GBARplus_{istate,jstate}(dd,dd);
    end
end
[TT,RR,gsig,theta_hat,retcode,options] = msre_solve(...
    AAminus_,AA0_,GBARplus_,BB_,...
    CC_,sparam,is_unique_steady_state,Q,T0,shock_horizon,options);

m_x=[];m_e=[];m_sig=[];
if ~retcode
    if any(stat_cols)
        m_sig=repmat({zeros(endo_nbr,1)},1,nstates);
        m_x=repmat({zeros(endo_nbr)},1,nstates);
        m_e=repmat({zeros(size(RR{1}))},1,nstates);
        solve_order=[orig_order(stat_cols),orig_order(dyn_cols)];
        stat_=1:static_endo_nbr;
        expect_order=size(RR{1},3);
        for istate=1:nstates
            iR11=AA0{istate}(stat_,stat_)\eye(static_endo_nbr);
            Di=AA0{istate}(stat_,dd);
            for jstate=1:nstates
                Di=Di+GBARplus{istate,jstate}(stat_,dd)*TT{jstate};
            end
            % risk
            m_sig{istate}(dd)=gsig{istate};
            tmp=Di*gsig{istate};
            for jstate=1:nstates
                tmp=tmp+GBARplus{istate,jstate}(stat_,dd)*gsig{jstate};
            end
            m_sig{istate}(stat_)=tmp;
            m_sig{istate}(solve_order)=m_sig{istate};
            
            % autoregressive part
            m_x{istate}(dd,dd)=TT{istate};
            m_x{istate}(stat_,dd)=-iR11*(AAminus{istate}(stat_,dd)+Di*TT{istate});
            m_x{istate}(solve_order,solve_order)=m_x{istate};
            
            % shock impact part
            m_e{istate}(dd,:,1)=RR{istate}(:,:,1);
            m_e{istate}(stat_,:,1)=-iR11*(BB{istate}(stat_,:)+Di*RR{istate}(:,:,1));
            for l=2:expect_order
                mei=Di*RR{istate}(:,:,l);
                for jstate=1:nstates
                    mei=mei+GBARplus{istate,jstate}(stat_,dd)*RR{jstate}(:,:,l-1);
                end
                m_e{istate}(stat_,:,l)=mei;
            end
            m_e{istate}(solve_order,:,:)=m_e{istate};
        end
    else
        m_x=TT;
        m_e=RR;
        m_sig=gsig;
    end
end


    function [AA0,AAminus,GBARplus,BB,CC,q]=msre_computational_savings()
        % this function separates static variables from dynamic ones. It places all
        % the static variables at the beginning and the dynamic ones following
        % after.
        % stat_cols is a logical vector indicating the original position of the
        % static variables.
        % the matrices are transformed in such a way that the dynamic equations
        % appear at the bottom and can be solved independently. Once the solution
        % for the dynamic variables is found, the solution for the static variables
        % can be computed as a function of the dynamic variables.
        
        AA0=A0;
        AAminus=Aminus;
        GBARplus=Gplus;
        BB=B;
        CC=C;
        q=cell(1,nstates);
        for st=1:nstates
            if any(stat_cols)
                [AA0{st},AAminus{st}]=...
                    re_order(A0{st},Aminus{st});
                [q{st},r]=qr(AA0{st}); %#ok<NASGU>
                AA0{st}=q{st}'*AA0{st};
                AAminus{st}=q{st}'*AAminus{st};
                BB{st}=q{st}'*BB{st};
                for slead=1:nstates
                    CC{st,slead}=q{st}'*CC{st,slead};
                    % Aplus still needs re-ordering before the application of the Q
                    % scheme
                    tmp=[GBARplus{st,slead}(:,stat_cols),GBARplus{st,slead}(:,~stat_cols)];
                    GBARplus{st,slead}=q{st}'*tmp;
                end
            else
                q{st}=1;
            end
        end
        function varargout=re_order(varargin)
            varargout=varargin;
            for ii=1:numel(varargin)
                varargout{ii}=[varargin{ii}(:,stat_cols),varargin{ii}(:,~stat_cols)];
            end
        end
    end
end