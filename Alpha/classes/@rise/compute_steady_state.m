function [obj,retcode]=compute_steady_state(obj,varargin)

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
    obj=struct('steady_state_file','');
    return
end

obj=set_options(obj,varargin{:});

steady_state_file=obj.options.steady_state_file;

endo_nbr=obj.endogenous.number(2);
% steady state functions
resid_func=@static_or_bgp_residuals;
last_item=endo_nbr;
func_ss=obj.func_handles.static;
func_jac=obj.model_derivatives.StaticEndogenous{1};
% func_jac=obj.func_handles.static_model_derivatives;
if ~isempty(obj.is_stationary_model)&& ~obj.is_stationary_model
    last_item=2*endo_nbr;
    func_ss=obj.func_handles.balanced_growth;
    func_jac=obj.model_derivatives.Static_BGP_Endogenous{1};
end

number_of_regimes=obj.markov_chains.regimes_number;
ss_and_bgp_start_vals=zeros(2*endo_nbr,number_of_regimes);

obj=compute_definitions(obj);
def=obj.solution.definitions;

x_ss=zeros(sum(obj.exogenous.number),1);

retcode=0;
if ~isempty(steady_state_file)
    if ischar(steady_state_file)
        steady_state_file=str2func(steady_state_file);
    end
    % this alternative is more general in the sense that it allows you to
    % modify the parameters, e.g. so that the model replicates the things
    % you want.
    nout=nargout(steady_state_file);
    if isempty(obj.steady_state_file_2_model_communication)
        % false means do not solve the steady state
        var_names=steady_state_file(obj,false);
        original_endo_ids=locate_variables(obj.orig_endo_names_current,var_names,true);
        located_vars_ids=find(~isnan(original_endo_ids));
        original_endo_ids=original_endo_ids(located_vars_ids);
        obj.steady_state_file_2_model_communication=...
            struct('original_endo_ids',original_endo_ids,...
            'located_vars_ids',located_vars_ids);
        % variables whose names are found in var_names get a steady state. this
        % includes the auxiliary variables, which are given in their
        % current-period form in the cell string orig_endo_names_current
        %         % take care of the variables that were logged in the model file
        %         %--------------------------------------------------------------
        %         logVars=obj.endogenous.name(obj.endogenous.is_log_var);
        %         if ~isempty(logVars)
        %             % inflate the list to include the auxiliary log vars
        %             log_var_ids=[];
        %             for ivar=1:numel(logVars)
        %                 bingo=find(strcmp(logVars{ivar},obj.orig_endo_names_current));
        %                 log_var_ids=[log_var_ids;bingo(:)]; %#ok<AGROW>
        %             end %locate_variables(obj.orig_endo_names_current,var_names,true);
        %             logVars=obj.orig_endo_names_current(log_var_ids);
        %             % now just as above, find the variables in the returned names
        %             ss_and_bgp_start_vals(logVars,:)=log(ss_and_bgp_start_vals(logVars,:));
        %         end
        
    end
    % now we actually solve the steady state. The second input argument is
    % set to true
    if nout==4
        [ss,obj,retcode,is_imposed]=steady_state_file(obj,true);
        obj.is_imposed_steady_state=is_imposed;
    else
        [ss,obj,retcode]=steady_state_file(obj,true);
    end
    if retcode
        retcode=1; % flag on the steady state
        return
    end
    % if the number of columns of ss does not correspond to the number of
    % regimes, replicate the columns
    deficient_cols=number_of_regimes-size(ss,2);
    if deficient_cols<0
        % In principle, this will never happen in a lifetime
        error([mfilename,':: steady state files returns more steady state columns than the number of regimes'])
    elseif deficient_cols>0
        ss=[ss,ss(:,ones(1,deficient_cols))];
    end
    % now put the content of ss into ss_and_bgp_start_vals
    original_endo_ids=obj.steady_state_file_2_model_communication.original_endo_ids;
    located_vars_ids=obj.steady_state_file_2_model_communication.located_vars_ids;
    ss_and_bgp_start_vals(located_vars_ids,:)=ss(original_endo_ids,:);
    % exogenous auxiliary variables have steady state zero and will
    % correspond to the nan elements in original_endo_ids
    
elseif obj.options.use_steady_state_model
    ssfunc=obj.func_handles.steady_state_model;
    if ~isempty(ssfunc)
        optimopt=optimset('display','none','maxiter',400);
        [ss,retcode]=steady_state_evaluation(optimopt,ssfunc);
        if retcode
            retcode=1; % flag on the steady state
            if obj(1).options.debug
                decipher_error(retcode)
            end
            return
        end
        if obj.is_optimal_policy_model
            if isempty(obj.steady_state_file_2_model_communication)
                % if this is to be used in this case, perhaps it should
                % change name too.
                original_endo_ids=find(~obj.endogenous.is_lagrange_multiplier);
                obj.steady_state_file_2_model_communication=...
                    struct('original_endo_ids',original_endo_ids,...
                    'located_vars_ids',[]);
            end
            tmp=zeros(endo_nbr,obj.markov_chains.regimes_number);
            tmp(obj.steady_state_file_2_model_communication.original_endo_ids,:)=ss;
            ss=tmp; clear tmp;
        end
        % now put the content of ss into ss_and_bgp_start_vals
        ss_and_bgp_start_vals(1:endo_nbr,:)=ss;
        % exogenous auxiliary variables have steady state zero
    end
end

% Now the parameters below may have been modified in the steady state file.
% for instance, nan may have been removed.
pp=obj.parameter_values;

% compute the unique steady state based on the ergodic distribution
%------------------------------------------------------------------
if obj.is_unique_steady_state
    [PAI00,retcode]=initial_markov_distribution(obj.solution.Q,1);
    if retcode
        if obj(1).options.debug
            decipher_error(retcode)
        end
        return
    end
    pp_i=0;
    def_i=[];
    if ~def{1}
        def_i=0;
    end
    for ii=1:number_of_regimes
        pp_i=pp_i+PAI00(ii)*pp(:,ii);
        if ~isempty(def{ii})
            def_i=def_i+PAI00(ii)*def{ii};
        end
    end
end
% check that the steady states computed above are actually the steady
% states. If not, then compute the steady state(s)
if ~obj.is_optimal_policy_model && ~obj.is_imposed_steady_state
    for ii=number_of_regimes:-1:1 % backward to optimize speed
        % if the initial guess solves the steady state then proceed. Else
        % try and improve the initial guess through fsolve.
        if ~obj.is_unique_steady_state
            def_i=def{ii};
            pp_i=pp(:,ii);
        end
        if isempty(obj.is_stationary_model)
            % determine stationarity
            [is_stationary,ss_and_bgp_start_vals(:,ii),retcode]=...
                determine_stationarity_status(ss_and_bgp_start_vals(:,ii),pp_i,def_i);
            if ~retcode
                obj.is_stationary_model=is_stationary;
                % based on that, set the functions...
                if ~obj.is_stationary_model
                    last_item=2*endo_nbr;
                end
            end
        else
            % solve the steady state
            % input list is always 'y'  'x'  'ss'  'param'  'def'  's0'  's1'
            [ss_and_bgp_start_vals(1:last_item,ii),retcode]=...
                solve_steady_state(ss_and_bgp_start_vals(1:last_item,ii),x_ss,[],pp_i,def_i,[],[],...
                resid_func,obj.is_linear_model,obj.options.optimset);
        end
        if retcode
            break
        end
        if obj.is_unique_steady_state
            ss_and_bgp_start_vals=ss_and_bgp_start_vals(:,ii*ones(1,number_of_regimes));
            break
        end
    end
end

for ii=1:number_of_regimes
    obj.solution.ss{ii}=sparse(ss_and_bgp_start_vals(1:endo_nbr,ii));
    obj.solution.bgp{ii}=sparse(ss_and_bgp_start_vals(endo_nbr+1:end,ii));
end

    function [is_stationary,ys,retcode]=determine_stationarity_status(ss_and_bgp_i,pp_i,def_i)
        
        is_stationary=[];
        % try stationarity
        
        ys=ss_and_bgp_i; % input list is always 'y'  'x'  'ss'  'param'  'def'  's0'  's1'
        [ys(1:endo_nbr),retcode]=solve_steady_state(ss_and_bgp_i(1:endo_nbr),x_ss,[],pp_i,def_i,[],[],...
            @static_or_bgp_residuals,obj.is_linear_model,obj.options.optimset);
        
        if retcode
            % try nonstationarity
            func_ss=obj.func_handles.balanced_growth;
            func_jac=obj.model_derivatives.Static_BGP_Endogenous{1};
            [ys,retcode]=solve_steady_state(ss_and_bgp_i,x_ss,[],pp_i,def_i,[],[],...
                @static_or_bgp_residuals,obj.is_linear_model,obj.options.optimset);
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

    function [lhs,Jac,retcode]=static_or_bgp_residuals(y_,x_,~,param_,def_,~,~)
        % if the structural model is defined in terms of steady states, then the
        % the word ss will appear in the static model. But then, those values are
        % the same as those found in y.
        % In some cases, the model is nonstationary so that only the BGP is
        % solvable in this case, the y vector also contains the balanced-growth
        % values. copying the whole vector to ss does not harm since only the
        % relevant elements in the upper part of y will be picked up.
        
        retcode=0;
        
        ss_=y_;
        lhs=online_function_evaluator(func_ss,y_,x_,ss_,param_,def_,[],[]);
        if nargout>1
            Jac=online_function_evaluator(func_jac,y_,x_,ss_,param_,def_,[],[]);
        end
    end

    function [ss,retcode]=steady_state_evaluation(options,ssfunc) %#ok<INUSL>
        
        options=optimset('display','none','maxiter',400,'tolfun',1e-6); %#ok<NASGU>
        
        retcode=0;
        ss=[];
        for ireg=1:number_of_regimes
            % in ssfunc, the definitions have been substituted already in
            % load_functions and so all is needed here is the parameters
            [y,obj.parameter_values(:,ireg)]=online_function_evaluator(...
                ssfunc,[],[],[],obj.parameter_values(:,ireg),[],[],[]); % y, x, ss, param, def, s0, s1
            retcode=~all(isfinite(y));
            if retcode
                return
            end
            if ireg==1
                ss=nan(numel(y),number_of_regimes);
            end
            ss(:,ireg)=y; %#ok<AGROW>
            if obj.is_unique_steady_state
                ss=ss(:,ireg*ones(1,number_of_regimes));
                break
            end
        end
    end
end
