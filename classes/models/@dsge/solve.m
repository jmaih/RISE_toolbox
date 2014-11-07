function [obj,retcode,structural_matrices]=solve(obj,varargin)
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
% obj=solve(obj,'solve_shock_horizon',struct('shock1',2,'shock3',4))
% obj=solve(obj,'solve_shock_horizon',5)
%
% See also:

% structrual_matrices not stored into the object in case we need to "get"
% it? in that case we enter the solving operations without solving and
% there is a flag to solve or not.

% I should formally check that all the variables that enter the calculation
% of endogenous probabilities have a unique steady state. I can add a flag
% into the endogenous that says whether a variables enters the endogenous
% probabilities or not. Instead of always loading the first column of the
% steady state when evaluating endogenous probabilities as I think I do in
% the beta version of RISE.

% anticipation_order=1
% solve_order = 1
% resolve if any of those has changed:
% - params,
% - solve_order: if a higher order has previously been solved, no need to >
% - anticipation_order, if a higher order has previously been solved, no need to >
% - solution algorithm,
% - etc...

% Please always respect the order of the inputs: y, x, ss, param, def

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=utils.miscellaneous.mergestructures(optimal_policy_solver_h(obj),...
        dsge_solver_h(obj),...
        struct('solver',[],...
        'solve_order',1,...
        'check_stability',true,...
        'solve_derivatives_type','symbolic',...%['symbolic','numerical','automatic']
        'solve_accelerate',false,...
        'solve_use_disc',false,...
        'solve_shock_horizon',[],...
        'vectorized_code_threshold',150));%
    return
end

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end
nobj=numel(obj);
if nobj>1
    retcode=nan(1,nobj);
    structural_matrices=cell(1,nobj);
    nout=nargout;
    for iobj=1:nobj
        if nout<3
            [obj(iobj),retcode(iobj)]=solve(obj(iobj),varargin{:});
        else
            [obj(iobj),retcode(iobj),structural_matrices{iobj}]=solve(obj(iobj),varargin{:});
        end
    end
    return
end

% initialize elements that will be used by different sub-functions
%-----------------------------------------------------------------
def=[];
ys=[];
ss=[];
xss=zeros(sum(obj.exogenous.number),1);
lead_lag_incidence=obj.lead_lag_incidence.after_solve;
endo_nbr=obj.endogenous.number(end);
the_leads=find(lead_lag_incidence(:,1)>0);
the_current=(1:endo_nbr)';
the_lags=find(lead_lag_incidence(:,3)>0);
indices=[the_leads;the_current;the_lags];
structural_matrices=[];

if obj.is_optimal_policy_model % || obj.estimation_under_way
    obj.options.solve_order=1;
end
solve_order=obj.options.solve_order;
params=obj.parameter_values;
resolve_it=check_whether_to_resolve();
h=obj.markov_chains.small_markov_chain_info.regimes_number;

% solve zeroth order or steady state ... and measurement errors
%--------------------------------------------------------------
retcode=solve_zeroth_order();
if solve_order>0 && ~retcode && resolve_it
    % the parameters may have changed during the computation of the zeroth
    % order or steady state. That is why they are reloaded here
    params=obj.parameter_values;
    
    % evaluate all the derivatives up to the desired order: analytical
    % derivatives can be evaluated sequentially, but for algorithmic
    % derivatives and for numerical derivatives, it is more economical to
    % do them at once, instead of calling the same functions several times
    %---------------------------------------------------------------------
    retcode=evaluate_all_derivatives();
    if ~retcode
        if obj.is_optimal_policy_model
            [T,eigval,zmultcols,retcode]=optimal_policy_solver_h(obj,structural_matrices);
            % expand solution to account for loose commitment
            %------------------------------------------------
            loose_com_col=find(strcmp(obj.markov_chains.chain_names,parser.loose_commit));
            if ~retcode && ~isempty(loose_com_col)
                big_regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));
                small_regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));
                bigh=obj.markov_chains.regimes_number;
                loose_com_regimes=big_regimes(:,loose_com_col);
                big_regimes(:,loose_com_col)=[];
                Tz=cell(1,bigh);
                % including the definitions, the steady state, etc.
                %--------------------------------------------------
                ss_=Tz;
                def_=Tz;
                bgp_=Tz;
                for ireg=1:bigh
                    bingo=all(bsxfun(@minus,big_regimes(ireg,:),small_regimes)==0,2);
                    Tsol=T.Tz{bingo};
                    if loose_com_regimes(ireg)==2 && ... % discretion: set multipliers to zero
                            ~obj.options.lc_reconvexify; % under reconvexification, we do not zero the multipliers
                        Tsol(:,zmultcols)=0;
                    end
                    Tz{ireg}=Tsol;
                    def_{ireg}=obj.solution.definitions{bingo};
                    ss_{ireg}=obj.solution.ss{bingo};
                    bgp_{ireg}=obj.solution.bgp{bingo};
                end
                obj.solution.definitions=def_;
                obj.solution.ss=ss_;
                obj.solution.bgp=bgp_;
                T.Tz=Tz;
            end
            % set the name of the solver
            %---------------------------
            if isempty(obj.options.solver)
                obj.options.solver='loose_commitment';
            end
        else
            [T,eigval,retcode,obj]=dsge_solver_h(obj,structural_matrices);
            % options may have changed and so we re-collect obj
        end
        inv_order_var=obj.inv_order_var.after_solve;
        if ~retcode
            T_fields=fieldnames(T);
            nregs=numel(T.Tz);
            for ifield_=1:numel(T_fields)
                fname=T_fields{ifield_};
                % re-order the rows right here, right now
                %----------------------------------------
                for ireg=1:nregs
                    T.(fname){ireg}=T.(fname){ireg}(inv_order_var,:);
                end
                obj.solution.(fname)=T.(fname);
            end
            obj.solution.eigval=eigval;
            
            obj=set_z_eplus_horizon(obj);
            obj.solution.user_resids=structural_matrices.user_resids;
        end
    end
end
if retcode && obj.options.debug
    utils.error.decipher(retcode)
end

    function retcode=evaluate_all_derivatives()
        % evaluate first-order derivatives
        %---------------------------------
        symbolic_type=strcmp(obj.options.solve_derivatives_type,'symbolic');
        xxx=repmat('v',1,solve_order);
        % evaluate higher-order derivatives
        %----------------------------------
        retcode=0;
        G01=cell(1,solve_order);
        for s1=1:h
            sparam=params(:,s1);
            for s0=1:h
                if ~retcode
                    % Note: G(s0,s1) =: ps0(s0,s1)*F(s0)
                    if symbolic_type
                        [G01{1:solve_order}]=utils.code.evaluate_functions(obj.routines.probs_times_dynamic_derivatives,...
                            ys(:,s0),xss,ss(:,s0),params(:,s0),sparam,def{s0},s0,s1);
                    else
                        G01=utils.code.evaluate_automatic_derivatives(...
                            obj.routines.symbolic.probs_times_dynamic,solve_order,...
                            ys(:,s0),xss,ss(:,s0),params(:,s0),sparam,def{s0},s0,s1);
                    end
                    if utils.error.valid(G01)
                        % use the derivatives Gi to build dv, dvv, dvvv, ...
                        %---------------------------------------------------
                        for io=1:solve_order
                            if ~symbolic_type
                                [nr,nc]=size(G01{io});
                                G01{io}=sparse(reshape(G01{io},nr,nc));
                            end
                            structural_matrices.(['d',xxx(1:io)]){s0,s1}=G01{io};%
                        end
                    else
                        retcode=2; % nans in jacobian
                    end
                end
            end
        end
        % Compute planner information first
        %----------------------------------
        if obj.is_optimal_policy_model || obj.is_optimal_simple_rule_model
            planner=struct(...
                'objective',{cell(1,h)},...
                'commitment',{cell(1,h)},...
                'discount',{cell(1,h)},...
                'weights',{cell(1,h)}...
                );
            orig_endo_nbr=obj.endogenous.number(1);
            for s0=1:h
                if ~retcode
                    lcd=utils.code.evaluate_functions(...
                        obj.routines.planner_loss_commitment_discount,...
                        ss(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s0);
                    if ~utils.error.valid(lcd)
                        retcode=6;
                    end
                    if ~retcode
                        % pick the hessian directly
                        order__=2;
                        if symbolic_type
                            [~,weights]=utils.code.evaluate_functions(...
                                obj.routines.planner_objective_derivatives,...
                                ss(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s0);
                        else
                            s1=s0;
                            D01=utils.code.evaluate_automatic_derivatives(...
                                obj.routines.symbolic.planner_objective,order__,...
                                ys(:,s0),xss,ss(:,s0),params(:,s0),[],def{s0},s0,s1);
                            weights=D01{order__};
                        end
                        if ~utils.error.valid(weights)
                            retcode=6;
                        end
                        if ~retcode
                            planner.objective{s0}=lcd(1);
                            planner.commitment{s0}=lcd(2);
                            planner.discount{s0}=lcd(3);
                            planner.weights{s0}=sparse(reshape(weights,orig_endo_nbr,orig_endo_nbr));
                        end
                    end
                end
            end
            structural_matrices.planner=planner;
            % change the ordering of the weight for the user. OSR will use
            % the structural matrices, which are ordered according to
            % order_var
            iov=obj.inv_order_var.before_solve;
            for s0=1:h
                planner.weights{s0}=planner.weights{s0}(iov,iov);
            end
            obj.solution.planner=planner; clear planner
        end
    end
    function resolve_it=check_whether_to_resolve()
        % forward expansion order
        % if order==0, it means there is no ancitipation. It does not mean
        % that there no shocks. If order==0, this will later on help set
        % the conditional information to nan. But for the purpose of
        % expanding the model, the expansion order is at least 1.
        % NORMALLY ALL OPTIONS SHOULD BE CHECKED. FOR INSTANCE ONE MAY WANT
        % TO SEE THE EFFECT OF CHANGING THE TOLERANCE LEVEL ON THE SOLUTION
        % AND CURRENTLY, THIS IS NOT TAKEN CARE OF. AT SOME LEVEL, ALL THE
        % OPTIONS COULD BE REPRESENTED IN A VECTOR, SO THAT IF AN ELEMENT
        % IS DIFFERENT, WHATEVER IT IS, THEN WE RESOLVE. BUT THIS IS NOT AN
        % URGENT MATTER.
        resolve_it=isempty(obj.current_solution_state)||...
            obj.estimation_under_way||...
            strcmp(obj.options.solve_initialization,'random')||...
            ~isfield(obj.solution,'Tz');
        horizon=max(obj.exogenous.shock_horizon);
        load_ssfuncs=isempty(obj.current_solution_state);
        if resolve_it
            obj.current_solution_state=fill_solve_state_info();
        else
            new_state=fill_solve_state_info();
            load_ssfuncs=load_ssfuncs||...
                ~strcmp(new_state.derivatives_type,obj.current_solution_state.derivatives_type)||...
                ~strcmp(new_state.solve_use_disc,obj.current_solution_state.solve_use_disc);
            resolve_it= ~isequal(new_state,obj.current_solution_state);
            if resolve_it
                obj.current_solution_state=new_state;
            end
        end
        % steady state functions (just for output)
        %-----------------------------------------
        if load_ssfuncs
            obj.steady_state_funcs=recreate_steady_state_functions(obj);
        end
        
        function new_state=fill_solve_state_info()
            new_state=struct('horizon',horizon,...
                'solution_algo',obj.options.solver,...
                'solve_order',solve_order,...
                'params',params,...
                'derivatives_type',obj.options.solve_derivatives_type,...
                'solve_use_disc',obj.options.solve_use_disc...
                );
        end
    end

    function retcode=solve_zeroth_order()
        retcode=0;
        if resolve_it
            obj.solution=struct();%dsge.initialize_solution_or_structure('solution',h);
            obj.solution.H=cell(1,h);
            [obj,structural_matrices,retcode]=compute_steady_state(obj);
            if ~retcode
                % measurement errors
                %-------------------
                if ~retcode
                    measure_flag=~isempty(obj.measurement_errors_restrictions);
                    if measure_flag
                        Restrictions=obj.measurement_errors_restrictions;
                        tmp=zeros(obj.observables.number(1),1);
                        for s0=1:h
                            tmp(Restrictions(:,1))=params(Restrictions(:,2),s0).^2;
                            obj.solution.H{s0}=diag(tmp);
                        end
                    end
                end
            end
        end
        % spit out the forms to be used for the computation of derivatives
        %-----------------------------------------------------------------
        if ~retcode
            ys=zeros(numel(indices),h);
            ss=zeros(endo_nbr,h);
            def=obj.solution.definitions;
            for s0=1:h
                for s1=1:h
                    if s0==1
                        ss(:,s1)=obj.solution.ss{s1};
                        ys(:,s1)=ss(indices,s1);
                        % y needs to be modified for the variables that have a
                        % nonzero balanced-growth path
                        if ~obj.is_stationary_model
                            if s1==1
                                [bgp_coefs,change_loc]=balanced_growth_path_powers();
                            end
                            ys(change_loc,s1)=ys(change_loc,s1).*obj.solution.bgp{s1}(indices(change_loc)).^bgp_coefs;
                        end
                    end
                end
            end
        end
        function [bgp_coefs,change_loc]=balanced_growth_path_powers()
            bgp_1=obj.solution.bgp{s1};
            bgp_vars=find(abs(bgp_1)>1e-3); % very low threshold...
            c_leads=zeros(size(the_leads)); c_leads(ismember(the_leads,bgp_vars))=1;
            c_current=zeros(size(the_current)); %c_current(ismember(the_current,bgp_vars))=0;
            c_lags=zeros(size(the_lags)); c_lags(ismember(the_lags,bgp_vars))=-1;
            cc=[c_leads(:);c_current(:);c_lags(:)];
            change_loc=cc~=0;
            bgp_coefs=cc(change_loc);
        end
    end
end

function ssfuncs=recreate_steady_state_functions(obj)
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