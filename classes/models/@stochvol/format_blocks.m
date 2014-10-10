function [endo_ids,exo_ids,params_ids,build_A,vectorize_A,build_OMG,vectorize_OMG]=format_blocks(obj)
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


% y-variables and y-shocks locations
%-----------------------------------
y_id=obj.endogenous_positions.y{2};
nobs=numel(y_id);
eps_y_id=obj.exogenous_positions.E{2};

% a-variables and a-shocks locations
%-----------------------------------
if obj.time_varying_parameters(1)
    % variables and shocks
    a_tmp=cell2mat(obj.endogenous_positions.a_{2}');
    a_state_id=a_tmp(:,1);
    a_matrix_id=a_tmp(:,3);
    eta_a_id=obj.exogenous_positions.ETA_a{2};
    % parameters
    a_tmp=cell2mat(obj.parameters_positions.a{2}');
    param_a_state_id=a_tmp(:,1);
    param_a_matrix_id=a_tmp(:,3);
    param_rho_a_state_id=obj.parameters_positions.rho_a{2};
    param_theta_a_state_id=obj.parameters_positions.theta_a{2};
else
    a_state_id=[];
    a_matrix_id=[];
    eta_a_id=[];
    param_a_state_id=[];
    param_a_matrix_id=[];
    param_rho_a_state_id=[];
    param_theta_a_state_id=[];
end

% sig-variables and sig-shocks locations
%-----------------------------------
if obj.time_varying_parameters(2)
    sig_id=obj.endogenous_positions.sig_{2};
    upsil_sig_id=obj.exogenous_positions.UPSIL_sig{2};
    % parameters
    param_sig_state_id=obj.parameters_positions.sig{2};
    param_rho_sig_state_id=obj.parameters_positions.rho_sig{2};
    param_theta_sig_state_id=obj.parameters_positions.theta_sig{2};
else
    sig_id=[];
    upsil_sig_id=[];
    param_sig_state_id=[];
    param_rho_sig_state_id=[];
    param_theta_sig_state_id=[];
end

% omg-variables and omg-shocks locations
%---------------------------------------
if obj.time_varying_parameters(3)
    omg_tmp=cell2mat(obj.endogenous_positions.omg_{2}');
    omg_state_id=omg_tmp(:,1);
    omg_matrix_id=omg_tmp(:,3);
    chi_omg_id=obj.exogenous_positions.CHI_omg{2};
    % parameters    
    omg_tmp=cell2mat(obj.parameters_positions.omg{2}');
    param_omg_state_id=omg_tmp(:,1);
    param_omg_matrix_id=omg_tmp(:,3);
    param_rho_omg_state_id=obj.parameters_positions.rho_omg{2};
    param_theta_omg_state_id=obj.parameters_positions.theta_omg{2};
else
    omg_state_id=[];
    omg_matrix_id=[];
    chi_omg_id=[];
    param_omg_state_id=[];
    param_omg_matrix_id=[];
    param_rho_omg_state_id=[];
    param_theta_omg_state_id=[];
end

% summing up
%-----------
endo_ids=struct('y_id',y_id,'a_state_id',a_state_id,'a_matrix_id',a_matrix_id,...
    'sig_id',sig_id,'omg_state_id',omg_state_id,'omg_matrix_id',omg_matrix_id);

exo_ids=struct('eps_y_id',eps_y_id,'eta_a_id',eta_a_id,'chi_omg_id',chi_omg_id,...
    'upsil_sig_id',upsil_sig_id);

params_ids=struct('a_state_id',param_a_state_id,'a_matrix_id',...
    param_a_matrix_id,'rho_a_state_id',param_rho_a_state_id,...
    'theta_a_state_id',param_theta_a_state_id,...
    'sig_state_id',param_sig_state_id,'rho_sig_state_id',...
    param_rho_sig_state_id,'theta_sig_state_id',param_theta_sig_state_id,...
    'omg_state_id',param_omg_state_id,'omg_matrix_id',param_omg_matrix_id,...
    'rho_omg_state_id',param_rho_omg_state_id,'theta_omg_state_id',...
    param_theta_omg_state_id);

% build up of matrices
%---------------------
build_A=@a_building;
vectorize_A=@a_vectorize;
build_OMG=@omg_building;
vectorize_OMG=@omg_vectorize;

    function ab=a_building(av)
        ab=zeros(nobs,nobs*obj.nlags+obj.nx);
        ab(a_matrix_id)=av;
    end
    function av=a_vectorize(ab)
        av=ab(a_matrix_id);
    end
    function ab=omg_building(av)
        ab=eye(nobs);
        ab(omg_matrix_id)=av;
    end
    function av=omg_vectorize(ab)
        av=ab(omg_matrix_id);
    end
end