function [x1,x1_linear]=simulation_engine(obj,...
    x0,... % initial conditions
    x0_linear,... % initial conditions for pruning
    z,... % deterministic variables
    istate,... % regime
    shocks)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

persistent simul_with_shocks is_observed exo_nbr endo_nbr ...
    A RHO_A THETA_A SIG RHO_SIG THETA_SIG OMG RHO_OMG THETA_OMG ...
    build_A vectorize_A build_OMG vectorize_OMG ...
    y_id eps_y_id a_state_id eta_a_id sig_id upsil_sig_id omg_state_id chi_omg_id
if nargin<6
    shocks=[];
end
if isempty(x0_linear)
    initialize_simulation();
end
if isempty(shocks)
    % the shocks include the deterministic variables
    %-----------------------------------------------
    shocks=generic_tools.set_exogenous_data(exo_nbr,is_observed,simul_with_shocks,z);
end
% initialize state vector for next period
%------------------------------------
x1=nan(endo_nbr,1);

% update Bt
%----------
At=vectorize_A(A{istate});
if obj.time_varying_parameters(1)
    % create new shocks eta for each time-varying parameter
    %------------------------------------------------------
    At_lag=x0(a_state_id,end);
    At=At+vectorize_A(RHO_A{istate}).*(At_lag-At);
    if size(shocks,2)
        eta_t=shocks(eta_a_id);
        At=At+vectorize_A(THETA_A{istate}).*eta_t;
    end
    % push into the state vector
    %---------------------------
    x1(a_state_id)=At;
end
% update SIGt
%--------------
SIGt=diag(SIG{istate});
if obj.time_varying_parameters(2)
    SIGt_lag=x0(sig_id,end);
    SIGt=SIGt.*(SIGt_lag./SIGt).^diag(RHO_SIG{istate});
    if size(shocks,2)
        upsil_sig_t=shocks(upsil_sig_id);
         SIGt=SIGt.*exp(diag(THETA_SIG{istate}).*upsil_sig_t);
    end
    % re-vectorize and push into the state vector
    %---------------------------------------------
    x1(sig_id)=SIGt;
end

% update OMGt
%------------
OMGt=vectorize_OMG(OMG{istate});
if obj.time_varying_parameters(3)
    OMGt_lag=x0(omg_state_id,end);
    OMGt=OMGt+vectorize_OMG(RHO_OMG{istate}).*(OMGt_lag-OMGt);
    if size(shocks,2)
        chi_omgt=shocks(chi_omg_id);
        OMGt=OMGt+vectorize_OMG(THETA_OMG{istate}).*chi_omgt;
    end
    % re-vectorize and push into the state vector
    %---------------------------------------------
    x1(omg_state_id)=OMGt;
end

% autoregressive part
%--------------------
x0_hat=x0(y_id,end:-1:1);
Xt=x0_hat(:);
% deterministic variables
%------------------------
Xt=[Xt;z];
if obj.constant
    Xt=[Xt;1];
end

% forecast y and store
%---------------------
y=build_A(At)*Xt;
if size(shocks,2)
    % create auxiliary v
    %-------------------
    e_t=shocks(eps_y_id);
    v=build_OMG(OMGt)*diag(SIGt)*e_t;
    y=y+v;
end
x1(y_id)=y;

% the full state vector
%----------------------
x1=[x0(:,2:end),x1];
x1_linear=x1;

    function initialize_simulation()
        simul_with_shocks=~obj.options.simul_no_shocks;
        is_observed=obj.exogenous.is_observed;
        exo_nbr=sum(obj.exogenous.number);
        endo_nbr=obj.endogenous.number(end);
        
        [A,RHO_A,THETA_A,SIG,RHO_SIG,THETA_SIG,OMG,RHO_OMG,THETA_OMG]=...
            stochvol.form_parameter_matrices(obj);
        
        [endo,exo,~,build_A,vectorize_A,build_OMG,vectorize_OMG]=...
            stochvol.format_blocks(obj);
        
        % endogenous locations
        %---------------------
        y_id=endo.y_id;
        a_state_id=endo.a_state_id;
%         a_matrix_id=endo.a_matrix_id;
        sig_id=endo.sig_id;
        omg_state_id=endo.omg_state_id;
%         omg_matrix_id=endo.omg_matrix_id;
        
        % shocks locations
        %-----------------
        eps_y_id=exo.eps_y_id;
        eta_a_id=exo.eta_a_id;
        chi_omg_id=exo.chi_omg_id;
        upsil_sig_id=exo.upsil_sig_id;
    end
end
