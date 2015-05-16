function [obj,ss_and_bgp_start_vals,retcode]=initial_steady_state(obj,varargin)
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


obj=set(obj,varargin{:});

number_of_regimes=obj.markov_chains.small_markov_chain_info.regimes_number;
endo_nbr=obj.endogenous.number(end);
ss_and_bgp_start_vals=zeros(2*endo_nbr,number_of_regimes);

[obj,retcode]=compute_definitions(obj);
if retcode
    return
end

steady_state_file=obj.options.steady_state_file;

if ~isempty(steady_state_file)
    if ischar(steady_state_file)
        steady_state_file=str2func(steady_state_file);
        obj.options.steady_state_file=steady_state_file;
    end
    % this alternative is more general in the sense that it allows you to
    % modify the parameters, e.g. so that the model replicates the things
    % you want.
    if isempty(obj.steady_state_file_2_model_communication)
        nout=nargout(steady_state_file);
        % false means do not solve the steady state
        var_names=steady_state_file(obj,false);
        original_endo_ids=locate_variables(obj.orig_endo_names_current,var_names,true);
        located_vars_ids=find(~isnan(original_endo_ids));
        original_endo_ids=original_endo_ids(located_vars_ids);
        obj.steady_state_file_2_model_communication=...
            struct('original_endo_ids',original_endo_ids,...
            'located_vars_ids',located_vars_ids,'nargout',nout);
    end
    % now we actually solve the steady state. The second input argument is
    % set to true
    if obj.steady_state_file_2_model_communication.nargout==4
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
    nx=sum(obj.exogenous.number);xss=zeros(nx,1);
    if ~isempty(ssfunc) && (isa(ssfunc,'function_handle')||...
            (isstruct(ssfunc)&&~isempty(ssfunc.code)))
        [ss,retcode,change_locs,pp_sstate]=steady_state_evaluation(ssfunc);
        if retcode
            retcode=1; % flag on the steady state
            if obj(1).options.debug
                utils.error.decipher(retcode)
            end
        else
            % now put the content of ss into ss_and_bgp_start_vals
            ss_and_bgp_start_vals(1:endo_nbr,:)=ss;
            % exogenous auxiliary variables have steady state zero
            
            % now replace the parameters
            %----------------------------
            if ~isempty(change_locs)
                obj.parameter_values(change_locs,:)=pp_sstate(change_locs,ones(1,number_of_regimes));
            end
        end
    end
end

if retcode && obj.options.debug
    utils.error.decipher(retcode)
end

    function [ss,retcode,change_locs,pp_sstate]=steady_state_evaluation(ssfunc)
        ss=[];
        change_locs=[];
        [pp_sstate,def_sstate,retcode]=get_sstate_parameters();
        if ~retcode
            pp_update=pp_sstate;
            for ireg=1:number_of_regimes
                [y,pp_update(:,ireg)]=utils.code.evaluate_functions(...
                    ssfunc,[],xss,[],pp_sstate(:,ireg),[],def_sstate{ireg},[],[]); % y, x, ss, param, sparam, def, s0, s1
                retcode=~(all(isfinite(y)) && isreal(y));
                if retcode
                    return
                end
                if ireg==1
                    ss=nan(numel(y),number_of_regimes);
                end
                y=auxiliary_endo_sstate_evaluation(obj,y,xss,pp_sstate(:,ireg),def_sstate{ireg});

                ss(:,ireg)=y; %#ok<AGROW>
                if obj.is_unique_steady_state
                    ss=ss(:,ireg*ones(1,number_of_regimes));
                    break
                end
            end
            % check whether the parameters have been updated and if so push
            %--------------------------------------------------------------
            first_col=pp_sstate(:,1);
            first_col_update=pp_update(:,1);
            nanlocs=isnan(first_col);
            if any(nanlocs(:))
                nanlocs_update=isnan(first_col_update);
                ddd=nanlocs-nanlocs_update;
                if any(ddd(:))
                    change_locs=nanlocs & ~nanlocs_update;
                    % the nans are not in the same place: the parameters have changed
                end
            end
            pp_sstate=pp_update;
        end
        
        function [pp_sstate,def_sstate,retcode]=get_sstate_parameters()
            obj=derive_auxiliary_parameters(obj);
            
            pp=obj.parameter_values;
            def=obj.solution.definitions;
            
            % compute the unique steady state based on the ergodic distribution
            %------------------------------------------------------------------
            def_sstate=def;
            pp_sstate=pp;
            retcode=0;
            if obj.is_unique_steady_state
                % Something here to improve upon: we compute the steady state based on
                % some information that we don't have yet, namely ys0. It seems both
                % should be computed simultaneously... in other words, nonlinearly.
                % This may destroy the ergodic distribution or even return some
                % errors... The problem does not occurr with constant probabilities or
                % when steady states are different
                if obj.is_endogenous_switching_model
                    error('please report this problem to junior.maih@gmail.com')
                end
                ys0_extended=ss_and_bgp_start_vals(:,1);
                [TransMat,retcode]=compute_steady_state_transition_matrix(...
                    obj.routines.transition_matrix,ys0_extended,pp(:,1),...
                    def{1},sum(obj.exogenous.number));
                if ~retcode
                    [pp_unique,def_unique,retcode]=...
                        dsge_tools.ergodic_parameters(TransMat.Qinit,def,pp);
                    % override the parameters and the definitions as they might be used
                    % for further processing in case the steady state is imposed
                    %------------------------------------------------------------------
                    pp_sstate=pp_unique(:,ones(1,number_of_regimes));
                    def_sstate=cellfun(@(x)def_unique,def_sstate,'uniformOutput',false);
                end
            end
        end
    end

end