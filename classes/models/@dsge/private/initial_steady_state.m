function [obj,ss_and_bgp_start_vals,ssfuncs,retcode]=initial_steady_state(obj,varargin)

obj=set(obj,varargin{:});

number_of_regimes=obj.markov_chains.small_markov_chain_info.regimes_number;
endo_nbr=obj.endogenous.number(end);
ss_and_bgp_start_vals=zeros(2*endo_nbr,number_of_regimes);

% steady state functions (just for output)
%-----------------------------------------
ssfuncs=recreate_steady_state_functions();

[obj,retcode]=compute_definitions(obj);
if retcode
    return
end

steady_state_file=obj.options.steady_state_file;

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
    else
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
    end
elseif obj.options.steady_state_use_steady_state_model
    ssfunc=obj.routines.steady_state_model;
    if ~isempty(ssfunc) %%%%%%%%%% && ~isempty(ssfunc.code)
        [ss,retcode]=steady_state_evaluation(ssfunc);
        if retcode
            retcode=1; % flag on the steady state
            if obj(1).options.debug
                utils.error.decipher(retcode)
            end
        else
            if obj.is_optimal_policy_model
                if isempty(obj.steady_state_file_2_model_communication)
                    % if this is to be used in this case, perhaps it should
                    % change name too.
                    original_endo_ids=find(~obj.endogenous.is_lagrange_multiplier);
                    obj.steady_state_file_2_model_communication=...
                        struct('original_endo_ids',original_endo_ids,...
                        'located_vars_ids',[]);
                end
                tmp=zeros(endo_nbr,number_of_regimes);
                tmp(obj.steady_state_file_2_model_communication.original_endo_ids,:)=ss;
                ss=tmp; clear tmp;
            end
            % now put the content of ss into ss_and_bgp_start_vals
            ss_and_bgp_start_vals(1:endo_nbr,:)=ss;
            % exogenous auxiliary variables have steady state zero
        end
    end
end

if retcode && obj.options.debug
    utils.error.decipher(retcode)
end

    function ssfuncs=recreate_steady_state_functions()
        % initialize this here
        ssfuncs=struct();
        
        symbolic_derivatives=strcmp(obj.options.solve_derivatives_type,'symbolic');
        ssfuncs.static=obj.routines.static;
        ssfuncs.static_bgp=obj.routines.static_bgp;
        if symbolic_derivatives
            ssfuncs.jac_static=@(varargin)utils.code.evaluate_functions(...
                obj.routines.static_derivatives,varargin{:});
            ssfuncs.jac_bgp=@(varargin)utils.code.evaluate_functions(...
                obj.routines.static_bgp_derivatives,varargin{:});
        else
            ssfuncs.jac_static=@(varargin)utils.code.compute_automatic_derivatives(...
                obj.routines.symbolic.static,1,...
                varargin{:});
            ssfuncs.jac_bgp=@(varargin)utils.code.compute_automatic_derivatives(...
                obj.routines.symbolic.static_bgp,1,...
                varargin{:});
        end
    end

    function [ss,retcode]=steady_state_evaluation(ssfunc) 
        
        options=optimset('display','none',...
            'maxiter',obj.options.fix_point_maxiter,...
            'tolfun',obj.options.fix_point_TolFun); %#ok<NASGU>
        
        retcode=0;
        ss=[];
        for ireg=1:number_of_regimes
            % in ssfunc, the definitions have been substituted already in
            % load_functions and so all is needed here is the parameters
            [y,obj.parameter_values(:,ireg)]=utils.code.evaluate_functions(...
                ssfunc,[],[],[],obj.parameter_values(:,ireg),[],[],[],[]); % y, x, ss, param, sparam, def, s0, s1
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