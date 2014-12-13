function [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also:


% In this file, there are two ways the steady state can be computed:
% case 1: a steady state function is passed as an argument
%---------------------------------------------------------
% the function should be of the form [params,ss,retcode,imposed]=myssfunc(params,flag)
% The inputs:
%           a-) params is the rise parameter object, which can be modified
%           inside the function and returned. In this way, some parameters
%           can be computed as functions of some pre-specified steady state
%           instead of the opposite.
%           b-) flag: when set to false, the function should return the
%           list of the variables for which the steady state is solved.
%           This gives to the user the flexibility to return the steady
%           state values in the order they want. When set to true, the
%           actual steady state should be computed.
% The outputs:
%           a-) params should be the rise parameter object, with or without
%           modifications
%           b-) ss should be the vector or matrix of the computed steady
%           states. If it is a matrix when there are multiple steady states
%           as could be the case in Markov switching. If you believe the
%           steady state is unique, you can still return a vector, which
%           will be multi-plicated if there are deficient columns.
%           c-) retcode should be 0 if there are no issues with the
%           computation of the steady state and 1 if the steady state
%           cannot be solved. Should you return a value different from 1,
%           no worries, RISE will correct the error for you.
% case 2: a steady state model is declared in the model file
%------------------------------------------------------------
% In the simplest case, the user does not have to worry about this
% function as long as the steady state model is correctly specified. The
% steady state model does not need to include all the variables of the
% model. The variables not included will automatically be assigned a
% starting value of 0 for the computation of the steady state. RISE will
% use the computed values in this block as starting values, should they not
% be the exact steady state.
% But this block could be rather complex as well, for instance with features like
% argzero, which would find the zero of sum function within the steady
% state model block.

% if definitions enter the solving of the steady state from an m-file
% written by the user, well, the user will just have to recompute them.

% one issue remains to be resolved, though, the fact that I change the
% configuration of the parameter object... The solution, perhaps would be
% to return to the old setting of keeping things the way they are in the
% original parameter object. But things will be slown down a smidge. In the
% long run, I should find a simple and efficient way of representing the
% parameters, which makes reading and writing easy.
% I also have to re-instate some old option steady_state_file


if isempty(obj)
    obj=struct('steady_state_file','',...
        'steady_state_use_steady_state_model',true);
    % ====== function handle that solves the steady state ======
    return
end

structural_matrices=struct();

[obj,ss_and_bgp_start_vals,retcode]=initial_steady_state(obj,varargin{:});

if retcode
    return
end

obj=derive_auxiliary_parameters(obj);

optimopt=obj.options.optimset;
optimopt.debug=obj.options.debug;
optimopt.is_unique=obj.is_unique_steady_state ;
optimopt.exo_nbr=sum(obj.exogenous.number);
optimopt.trans_mat_func=obj.routines.transition_matrix;
optimopt.is_linear_model=obj.is_linear_model;
x_ss=zeros(sum(obj.exogenous.number),1);
pp=obj.parameter_values;
def=obj.solution.definitions;

% don't look for this elsewhere. It has to be consistent with the above
%----------------------------------------------------------------------
number_of_regimes=size(ss_and_bgp_start_vals,2);

% initial steady state functions in case of stationarity
%-------------------------------------------------------
func_ss=obj.steady_state_funcs.static; % func_ss_static func_ss_bgp
func_jac=obj.steady_state_funcs.jac_static; % func_jac_static func_jac_bgp

endo_nbr=obj.endogenous.number(end);
if isempty(obj.is_stationary_model)|| obj.is_stationary_model
    last_item=endo_nbr;
else
    last_item=2*endo_nbr;
end
structural_matrices.user_resids=sparse(endo_nbr,number_of_regimes);

% compute the unique steady state based on the ergodic distribution
%------------------------------------------------------------------
def_sstate=def;
pp_sstate=pp;
if optimopt.is_unique
    % Something here to improve upon: we compute the steady state based on
    % some information that we don't have yet, namely ys0. It seems both
    % should be computed simultaneously... in other words, nonlinearly.
    % This may destroy the ergodic distribution or even return some
    % errors... The problem does not occurr with constant probabilities or
    % when steady states are different
    ys0=ss_and_bgp_start_vals(1:last_item,1);
    [TransMat,retcode]=compute_steady_state_transition_matrix(...
        optimopt.trans_mat_func,ys0,pp(:,1),def{1},optimopt.exo_nbr);
    if ~retcode
        [pp_unique,def_unique,retcode]=...
            dsge_tools.ergodic_parameters(TransMat.Qinit,def,pp);
        % override the parameters and the definitions as they might be used
        % for further processing in case the steady state is imposed
        %------------------------------------------------------------------
        pp_sstate=pp_unique(:,ones(1,number_of_regimes));
        def_sstate=cellfun(@(x)def_unique,def_sstate,'uniformOutput',false);
    end
    if retcode
        return
    end
end

if obj.is_optimal_policy_model||obj.is_imposed_steady_state
    if isempty(obj.is_stationary_model)
        obj.is_stationary_model=true;
        % this will have to change some day...
    end
else
    if isempty(obj.is_stationary_model)
        [is_stationary,ys,retcode]=determine_stationarity_status(ss_and_bgp_start_vals);
        if ~retcode
            obj.is_stationary_model=is_stationary;
            % based on that, set the functions...
            if ~obj.is_stationary_model
                last_item=2*endo_nbr;
            end
        end
    else
        [ys,retcode]=solve_steady_state(ss_and_bgp_start_vals(1:last_item,:),...
            def_sstate,pp_sstate,@ss_residuals,optimopt);
    end
    ss_and_bgp_start_vals(1:last_item,:)=ys(1:last_item,:);
end

if ~retcode
    if obj.is_imposed_steady_state
        % compute the residuals
        %----------------------
        user_resids=nan(endo_nbr,number_of_regimes);
        for ireg=1:number_of_regimes
            user_resids(:,ireg)=ss_residuals(ss_and_bgp_start_vals(1:last_item,ireg),...
                pp(:,ireg),def{ireg});
        end
        user_resids(abs(user_resids)<1e-9)=0;
        structural_matrices.user_resids=sparse(user_resids);
    end
    
    [obj.solution.transition_matrices,retcode]=...
        compute_steady_state_transition_matrix(obj.routines.transition_matrix,...
        ss_and_bgp_start_vals(:,1),pp(:,1),def{1},...
        sum(obj.exogenous.number));
    
    structural_matrices.transition_matrices=obj.solution.transition_matrices;
    
    % compute theta_hat depending on the choice of steady state
    %-----------------------------------------------------------
    
    theta_hat=cell(number_of_regimes);
    for s1=1:number_of_regimes
        for s0=1:number_of_regimes
            tmp=pp(:,s1)-pp_sstate(:,s0);
            theta_hat{s0,s1}=(1-obj.options.solve_disable_theta)*tmp(obj.steady_state_index.theta_plus);
        end
    end
    structural_matrices.theta_hat=theta_hat;
    
    ss_=ss_and_bgp_start_vals(1:endo_nbr,:);
    bgp_=ss_and_bgp_start_vals(endo_nbr+1:end,:);
    ss_tvp=ss_(obj.endogenous.is_affect_trans_probs,:);
    bad=any(abs(bsxfun(@minus,ss_tvp,ss_tvp(:,1)))>1e-9,2);
    if any(bad)
        bad_endo_vars=get(obj,'endo_list(affect_trans_probs)');
        bad_endo_vars=bad_endo_vars(bad);
        disp(bad_endo_vars)
        error(['The variables above affect the transition probabilities but ',...
            'do not have the same steady state in each regime'])
    end
    for ii=1:number_of_regimes
        obj.solution.ss{ii}=sparse(ss_(:,ii));
        obj.solution.bgp{ii}=sparse(bgp_(:,ii));
    end
end

    function [is_stationary,ys,retcode]=determine_stationarity_status(ss_and_bgp)
        
        is_stationary=[];
        % try stationarity
        
        ys=ss_and_bgp;
        [ys(1:endo_nbr,:),retcode]=solve_steady_state(ss_and_bgp(1:endo_nbr,:),...
            def_sstate,pp_sstate,@ss_residuals,optimopt);
        
        if retcode
            % try nonstationarity
            func_ss=obj.steady_state_funcs.static_bgp;
            func_jac=obj.steady_state_funcs.jac_bgp;
            
            [ys,retcode]=solve_steady_state(ss_and_bgp,...
                def_sstate,pp_sstate,@ss_residuals,optimopt);
            
            if ~retcode
                is_stationary=false;
                disp([mfilename,':: model is not mean-stationary but allows for a BALANCED GROWTH PATH'])
            else
                disp([mfilename,':: stationarity status could not be determined'])
            end
        else
            is_stationary=true;
            disp([mfilename,':: model found to be MEAN-stationary'])
        end
    end

    function [r,Jac,retcode]=ss_residuals(ss_i,pp_i,def_i)%func_ss,func_jac,x_ss,
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
        % More About
        % ------------
        %
        % Examples
        % ---------
        %
        % See also:
        
        
        retcode=0;
        Jac=[];
        
        y_=ss_i;
        % input list is always 'y'    'x'    'ss'    'param'    'sparam'    'def'    's0'    's1'
        r=utils.code.evaluate_functions(func_ss,y_,x_ss,ss_i,pp_i,[],def_i,[],[]);
        if nargout>1
            Jac=utils.code.evaluate_functions(func_jac,y_,x_ss,ss_i,pp_i,[],def_i,[],[]);
        end
        
        if retcode && obj(1).options.debug
            utils.error.decipher(retcode)
        end
    end
end
