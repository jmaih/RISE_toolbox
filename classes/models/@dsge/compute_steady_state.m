function [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)

% compute_steady_state - computes the steady state of a dsge model
%
% Syntax
% -------
% ::
%
%   [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model file
%
% - **varargin** []: usual optional arguments
%
% Outputs
% --------
%
% - **obj** [rise|dsge]: model file
%
% - **structural_matrices** [struct]: structure containing various
%   important elements for the solution of the system
%
% - **retcode** [numeric]: 0 if there was problem computing the steady
%   state.
%
% More About
% ------------
%
% - There are 2 cases to consider:
%   - The user does not provide any steady state equations: RISE will
%   attempt to solve the steady state using a vector of zeros as initial
%   guess. In this first step, RISE assumes the model is stationary. If
%   this steps fails, RISE assumes the process is nonstationary and
%   proceeds to solving for the balanced growth path.
%   - The user provide some equations for solving the steady state. This is
%   done in two ways:
%       1-) the steady_state_model block: the variables that do not appear
%           in the block will be initialized at zero. Some parameters can
%           also be computed inside the block. The user can define an
%           optimization to solve for a subset of steady state values
%           inside the block. The block has three attributes 
%           (i) imposed(default=false): RISE computes the solution at the
%           specified point without checking that the point solves for the
%           steady state
%           (ii) unique (default=false): RISE computes the steady state at
%           the ergodic distribution of the parameters. In case the
%           probabilities are endogenous, the ergodic distribution of the
%           parameters is itself a function of the steady state of the
%           variables.
%           (iii) initial_guess(default=false): If the calculation of the
%           steady state failes at the initial guess, RISE ignores the
%           equations in all subsequent iterations for finding the steady
%           state.
%       2-) the steady state file: The user writes a function which can be 
%           called in two possible ways: (i) [vnames,info]=ssfile(); In
%           this case the first output argument is the list of variables
%           for which the user computes the steady state; the second output
%           is a structure with fields unique, imposed and initial_guess
%           just as in the case of the steady state model. (ii) The other
%           call to the function is [y,newp,retcode]=ssfile(y,p,d,id,obj).
%           In this case, the first input (y) is the vector of steady
%           states, which is updated and returned as the first output. The
%           locations of the modifications are indicated by the fourth
%           input (id), which is computed based on the list of the
%           variables in vnames above. As for the other outputs, p is a
%           structure with parameters, d is a structure with definitions,
%           obj is the model object in case the user needs some further
%           information for computing the steady state. In case some
%           parameters are computed in the steady state file, they should
%           be returned in the structure "newp". The last output "retcode"
%           indicates whether no problem was encoutered in the computation
%           of the steady state (retcode=0) or the opposite (retcode =
%           any number different from 0).
%
% - Writing the steady state file in this ways makes it possible to use the
%   same whether there is regime switching or not.
%
% Examples
% ---------
%
% See also:

if isempty(obj)
    obj=struct('steady_state_file','',...
        'steady_state_use_steady_state_model',true,...
        'steady_state_solver','fsolve');
    % ====== function handle that solves the steady state ======
    return
end

structural_matrices=struct();

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end

number_of_regimes=obj.markov_chains.small_markov_chain_info.regimes_number;
endo_nbr=obj.endogenous.number;
exo_nbr=sum(obj.exogenous.number);
x=zeros(exo_nbr,1);
ss_bgp=zeros(2*endo_nbr,number_of_regimes);

[obj,retcode]=compute_definitions(obj);
if retcode
    return
end

optimopt=obj.options.optimset;
tol=optimopt.TolFun;
sqrt_tol=sqrt(tol);
debug=obj.options.debug;
arg_zero_solver=obj.options.steady_state_solver;

p=obj.parameter_values;
d=obj.solution.definitions;
pnames=obj.parameters.name;
dnames=obj.definitions.name;

[steady_state_file,ids,ssfile_solved,is_param_changed_in_ssfile]=...
    prepare_steady_state_program();
feasible_steady_state_file=~isempty(steady_state_file);

steady_state_model=obj.routines.steady_state_model;
is_param_changed_in_ssmodel=obj.is_param_changed_in_ssmodel;

feasible_steady_state_model=~feasible_steady_state_file &&...
    ~isempty(steady_state_model) && ...
    obj.options.steady_state_use_steady_state_model;

% after having parsed the ssfile and the ssmodel, we can do this
%----------------------------------------------------------------
is_initial_guess_steady_state=obj.is_initial_guess_steady_state;
errmsg=['With an initial guess for the steady state, ',...
        'the parameters cannot be modified in the steady state ',...
        'model or file'];
if is_initial_guess_steady_state && any(is_param_changed_in_ssmodel)
    error(errmsg)
end

if isempty(obj.is_stationary_model)
    if ~is_initial_guess_steady_state && ...
            (~isempty(steady_state_file)||~isempty(steady_state_model))
        obj.is_stationary_model=true;
    end
end

if isempty(obj.is_stationary_model)||obj.is_stationary_model
    static_model=obj.steady_state_funcs.static;
    static_model_jacobian=obj.steady_state_funcs.jac_static;
    y=ss_bgp(1:endo_nbr,:);
elseif ~isempty(obj.is_stationary_model) && ~obj.is_stationary_model
    static_model=obj.steady_state_funcs.static_bgp;
    static_model_jacobian=obj.steady_state_funcs.jac_bgp;
    y=ss_bgp;
end

[ss,r,retcode]=run_one_pass(y);

if ~retcode
    if isempty(obj.is_stationary_model)
        obj.is_stationary_model=true;
    end
else
    if isempty(obj.is_stationary_model)
        % stationarity failed above. Now investigate nonstationarity
        %------------------------------------------------------------
        static_model=obj.steady_state_funcs.static_bgp;
        static_model_jacobian=obj.steady_state_funcs.jac_bgp;
        y=ss_bgp;
        [ss,r,retcode]=run_one_pass(y);
        if ~retcode
            obj.is_stationary_model=false;
            % else stationarity status cannot be determined...
        end
    end
end

if ~retcode
    % update the parameters
    %-----------------------
    obj.parameter_values=p;
    
    % prepare output
    %-----------------
    structural_matrices.user_resids=r(1:endo_nbr,:);
    if obj.is_unique_steady_state
        [~,~,obj.solution.transition_matrices,retcode]=ergodic_parameters(ss(:,1));
    else
        [obj.solution.transition_matrices,retcode]=...
            compute_steady_state_transition_matrix(obj.routines.transition_matrix,...
            ss(:,1),p(:,1),d{1},sum(obj.exogenous.number));
    end
    structural_matrices.transition_matrices=obj.solution.transition_matrices;
    if obj.is_stationary_model
        ss_bgp(1:endo_nbr,:)=ss;
    else
        ss_bgp=ss;
    end
    ss_bgp=sparse(ss_bgp);
    ss_=ss_bgp(1:endo_nbr,:);
    bgp_=ss_bgp(endo_nbr+1:end,:);
    ss_tvp=ss_(obj.endogenous.is_affect_trans_probs,:);
    bad=any(abs(bsxfun(@minus,ss_tvp,ss_tvp(:,1)))>1e-9,2);
    if any(bad)
        bad_endo_vars=get(obj,'endo_list(affect_trans_probs)');
        bad_endo_vars=bad_endo_vars(bad);
        disp(bad_endo_vars)
        error(['The variables above affect the transition probabilities but ',...
            'do not have the same steady state in each regime'])
    end
    obj.solution.ss=mat2cell(ss_,endo_nbr,ones(1,number_of_regimes));
    obj.solution.bgp=mat2cell(bgp_,endo_nbr,ones(1,number_of_regimes));
end

%-------------------------------------------------------------------------%
% -------------------------- Nested functions --------------------------- %
%-------------------------------------------------------------------------%

    function flag=success(r)        
        flag=utils.error.valid(r) && solves_steady_state();
        function flag=solves_steady_state()
            flag=obj.is_imposed_steady_state||max(abs(r(:)))<sqrt_tol;
        end
    end

    function [ss,r,retcode]=run_one_pass(y)
        retcode=0;
        ss=y;
        r=y;
        nregs=obj.is_unique_steady_state*(1-number_of_regimes)+...
            number_of_regimes;
        for ireg=1:nregs
            p_ireg=p(:,ireg);
            [ss_1,p_ireg,retcode,r_1,Jac]=do_one_regime(y(:,ireg),p_ireg,...
                d{ireg});
            if retcode
                return
            else
                ss(:,ireg)=ss_1;
                r(:,ireg)=r_1;
                if success(r(:,ireg))
                    p(:,ireg)=p_ireg;
                else
                    if obj.is_linear_model
                        ss0=ss(:,ireg)-Jac\r(:,ireg);
                        [ss(:,ireg),p_ireg,retcode,r(:,ireg)]=do_one_regime(...
                            ss0,p_ireg,d{ireg});
                    else
                        vlist=~obj.auxiliary_variables.model;
                        if obj.is_initial_guess_steady_state
                            feasible_steady_state_file=false;
                            feasible_steady_state_model=false;
                        else
                            % instruments may come from ssmodel or ssfile
                            %---------------------------------------------
                            if feasible_steady_state_model
                                vlist=~obj.auxiliary_variables.ssmodel_solved & vlist;
                                %vlist=~obj.auxiliary_variables.sstate_solved;
                            elseif feasible_steady_state_file
                                vlist=~ssfile_solved & vlist;
                            end
                        end
                        missing=size(ss,1)-numel(vlist);
                        vlist=[vlist,true(1,missing)]; %#ok<AGROW>
                        if any(vlist)
                            % this can only be done if there are non-solved
                            % variables
                            [ss(:,ireg),r(:,ireg),retcode]=...
                                optimizer_over_vlist(vlist,ireg,ss(vlist,ireg));
                        else
                            % we simply return with the retcode
                        end
                    end
                    if ~retcode
                        if success(r(:,ireg))
                            if nregs<number_of_regimes
                                if feasible_steady_state_model
                                    p(is_param_changed_in_ssmodel,2:number_of_regimes)=p(is_param_changed_in_ssmodel,1)*ones(1,number_of_regimes-1);
                                elseif feasible_steady_state_file
                                    p(is_param_changed_in_ssfile,2:number_of_regimes)=p(is_param_changed_in_ssfile,1)*ones(1,number_of_regimes-1);
                                end
                            else
                                p(:,ireg)=p_ireg;
                            end
                        else
                            retcode=1;
                        end
                    end
                end
            end
            if retcode
                ss=[];r=[];
                return
            end
        end
        if number_of_regimes>1 && nregs<number_of_regimes
            ss=ss(:,ones(1,number_of_regimes));
            r=r(:,ones(1,number_of_regimes));
        end
        
        function [ss,resids,retcode]=optimizer_over_vlist(vlist,ireg,ssc0)
            if nargin<3
                ssc0=zeros(nv,1);
            end
            nv=sum(vlist);
            % check that we do not get problems starting at zero
            %---------------------------------------------------
            [~,~,~,rcode]=concentrated_residuals(ssc0);
            if rcode
                ssc0=rand(nv,1);
                if debug
                    warning('randomized initial conditions for sstate solver')
                end
            end
            if debug
                optimopt.Display='iter';
                %                 optimopt.Jacobian='off';
                %                 optimopt.MaxIter=100000;
                %                 optimopt.MaxFunEvals=10000;
            else
                optimopt.Display='none';
            end
            %--------------------------------------------------------------
            if strcmp(arg_zero_solver,'lsqnonlin')
                % call to lsqnonlin
                %------------------
                [ssc1,resnorm,residuals,exitflag]=lsqnonlin(...
                    @concentrated_residuals,ssc0,[],[],optimopt);  %#ok<ASGLU>
            elseif strcmp(arg_zero_solver,'fsolve')
                % call to fsolve
                %------------------
                [ssc1,fval,exitflag]=fsolve(@concentrated_residuals,ssc0,optimopt);
                resnorm=norm(fval);
            elseif strcmp(arg_zero_solver,'fminunc')
                % call to fsolve
                %------------------
                [ssc1,resnorm,exitflag]=fminunc(@minimizer,ssc0,optimopt);
            else
                [arg_zero_solver,vargs]=utils.code.user_function_to_rise_function(arg_zero_solver);
                [ssc1,resnorm,exitflag]=arg_zero_solver(@concentrated_residuals,ssc0,optimopt,vargs{:});
            end
            % N.B: I use sqrt(tol) below on purpose!!!
            %-----------------------------------------
            exitflag=utils.optim.exitflag(exitflag,ssc1,resnorm,sqrt_tol);
            %--------------------------------------------------------------
            retcode=1-(exitflag==1);
            [resids,Jac,ss]=concentrated_residuals(ssc1);
            function fval=minimizer(ssc)
                fval=concentrated_residuals(ssc);
                fval=norm(fval);
            end
            function [resids,Jac,ss1,rcode]=concentrated_residuals(ssc)
                sstmp=y(:,ireg);
                sstmp(vlist)=ssc;
                [ss1,p_ireg,rcode,resids,Jac]=do_one_regime(sstmp,p_ireg,...
                    d{ireg});
                if rcode||any(isnan(resids))||any(isinf(resids))
                    tmp=sqrt(obj.options.estim_penalty/(endo_nbr*nv));
                    Jac=tmp*ones(endo_nbr,nv);
                    resids=Jac(:,1);
                    % explicitly signal there is a problem
                    %--------------------------------------
                    rcode=inf; 
                else
                    Jac=Jac(:,vlist);
                end
            end
        end
        
        function [ss,pp,retcode,r,Jac]=do_one_regime(ys,pp,dd)
            retcode=0;
            ss=[];r=[];Jac=[];
            if obj.is_unique_steady_state
                [pp,dd,~,retcode]=ergodic_parameters(ys);
                if retcode,return,end
            end
            % compute first pass of steady state
            %------------------------------------
            if feasible_steady_state_file
                % here the parameters must be turned into a structure
                %-----------------------------------------------------
                [ss,newp,retcode]=steady_state_file(obj,...
                    ys,...
                    vector2struct(pnames,pp),...
                    vector2struct(dnames,dd),...
                    ids);
                
                if ~isempty(newp)
                    fnames=fieldnames(newp);
                    locs=locate_variables(fnames,obj.parameters.name);
                    if ~any(is_param_changed_in_ssfile)
                        is_param_changed_in_ssfile(locs)=true;
                    end
                    if is_initial_guess_steady_state && any(is_param_changed_in_ssfile)
                        error(errmsg)
                    end
                    for ipos=1:numel(locs)
                        pp(locs(ipos))=newp.(fnames{ipos});
                    end
                end
            elseif feasible_steady_state_model
                % here the parameters are a vector
                %----------------------------------
                % y, x, ss, param, def, s0, s1
                % with ys as input, vector ss will not be initialized!!!
                [ss,pp]=utils.code.evaluate_functions(steady_state_model,ys,x,[],pp,dd,[],[]);
            else
                ss=ys;
            end
            
            if retcode,return,end
                
            % compute auxiliary steady state
            %---------------------------------
            ss=auxiliary_endo_sstate_evaluation(obj,ss,x,pp,dd);
            
            % evaluate residuals
            %--------------------
            y_=ss;
            r=utils.code.evaluate_functions(static_model,y_,x,ss,pp,dd,[],[]); % func_ss
            Jac=utils.code.evaluate_functions(static_model_jacobian,y_,x,ss,pp,dd,[],[]);%func_jac
        end
    end

    function [steady_state_file,ids,ssfile_solved,is_param_changed_in_ssfile]=prepare_steady_state_program()
        is_param_changed_in_ssfile=false(1,sum(obj.parameters.number));
        if ~isempty(obj.options.steady_state_file) && ...
                ~isa(obj.options.steady_state_file,'function_handle')
            if ~ischar(obj.options.steady_state_file)
                error('steady state file must be a function handle or a char')
            end
            obj.options.steady_state_file=str2func(obj.options.steady_state_file);
        end
        steady_state_file=obj.options.steady_state_file;
        ids=[];
        ssfile_solved=[];
        if ~isempty(steady_state_file)
            if isempty(obj.steady_state_file_2_model_communication)
                % call the function with the object only
                %----------------------------------------
                [var_names,new_settings]=steady_state_file(obj);
                default_settings=struct('unique',false,'imposed',false,...
                    'initial_guess',true);
                ff=fieldnames(new_settings);
                for ifield=1:numel(ff)
                    if ~isfield(default_settings,ff{ifield})
                        error(['The expected fields of a steady state file are : ',...
                            '"unique","imposed" and "initial_guess"'])
                    end
                    default_settings.(ff{ifield})=new_settings.(ff{ifield});
                end
                obj.is_imposed_steady_state=default_settings.imposed;
                obj.is_unique_steady_state=default_settings.unique;
                obj.is_initial_guess_steady_state=default_settings.initial_guess;
                user_endo_ids=locate_variables(var_names,obj.endogenous.name);
                ssfile_solved=ismember(obj.endogenous.name,var_names);
                obj.steady_state_file_2_model_communication=...
                    struct('user_endo_ids',user_endo_ids,...
                    'ssfile_solved',ssfile_solved);
            end
            ids=obj.steady_state_file_2_model_communication.user_endo_ids;
            ssfile_solved=obj.steady_state_file_2_model_communication.ssfile_solved;
        end
    end

    function [pp_unique,def_unique,TransMat,retcode]=ergodic_parameters(y)
        [TransMat,retcode]=compute_steady_state_transition_matrix(...
            obj.routines.transition_matrix,y,p(:,1),...
            d{1},sum(obj.exogenous.number));
        if retcode
            pp_unique=[];
            def_unique=[];
        else
            [pp_unique,def_unique,retcode]=...
                dsge_tools.ergodic_parameters(TransMat.Qinit,d,p);
        end
    end
end

%-------------------------------------------------------------------------%
% ---------------------------- Sub-functions ---------------------------- %
%-------------------------------------------------------------------------%

function pnew=vector2struct(vnames,v)
pnew=struct();
for ip=1:numel(v)
    pnew.(vnames{ip})=v(ip);
end
end