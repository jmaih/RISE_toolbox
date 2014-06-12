function [obj,retcode,structural_matrices]=solve(obj,varargin)
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
                    if loose_com_regimes(ireg)==2 % discretion: set multipliers to zero
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
        else
            [T,eigval,retcode]=dsge_solver_h(obj,structural_matrices);
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
            if symbolic_type
                for io=1:solve_order
                    [Gi,retcode]=evaluate_structural_derivatives(obj.routines.probs_times_dynamic_derivatives(io));
                    if ~retcode
                        % use the derivatives Gi to build dv, dvv, dvvv, ...
                        %---------------------------------------------------
                        for s0=1:h
                            for s1=1:h
                                structural_matrices.(['d',xxx(1:io)]){s0,s1}=Gi{s0,s1};%
                            end
                        end
                    end
                end
            else
                retcode=0;
                sparam=params;
                for s0=1:h
                    for s1=1:h
                        if ~retcode
                            [D01,retcode]=utils.code.compute_automatic_derivatives(...
                                obj.routines.symbolic.probs_times_dynamic,solve_order,...
                                ys(:,s0),xss,ss(:,s0),params(:,s0),sparam,def{s0},s0,s1);
                            for io=1:solve_order
                                [nr,nc]=size(D01{io});
                                structural_matrices.(['d',xxx(1:io)]){s0,s1}=...
                                     sparse(reshape(D01{io},nr,nc));
                            end
                        end
                    end
                end
                clear sparam
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
                        if symbolic_type
                            % pick the hessian directly
                            weights=utils.code.evaluate_functions(...
                                obj.routines.planner_objective_derivatives(2),...
                                ss(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s0);
                        else
                            up_to_hessian=2;
                            s1=s0;
                            [D01,retcode]=utils.code.compute_automatic_derivatives(...
                                obj.routines.symbolic.planner_objective,up_to_hessian,...
                                ys(:,s0),xss,ss(:,s0),params(:,s0),[],def{s0},s0,s1);
                            weights=D01{2};
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
        if resolve_it
            obj.current_solution_state=fill_solve_state_info();
        else
            new_state=fill_solve_state_info();
            resolve_it= ~isequal(new_state,obj.current_solution_state);
            if resolve_it
                obj.current_solution_state=new_state;
            end
        end
        function new_state=fill_solve_state_info()
            new_state=struct('horizon',horizon,...
                'solution_algo',obj.options.solver,...
                'solve_order',solve_order,...
                'params',params,...
                'derivatives_type',obj.options.solve_derivatives_type...
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

    function [D,retcode]=evaluate_structural_derivatives(derivatives)
        retcode=0;
        D=cell(h);
        sparam=params;
        for s0=1:h
            for s1=1:h
                % Note: G(s0,s1) =: ps0(s0,s1)*F(s0)
                D{s0,s1}=utils.code.evaluate_functions(derivatives,...
                    ys(:,s0),xss,ss(:,s0),params(:,s0),sparam,def{s0},s0,s1);
                if ~utils.error.valid(D{s0,s1})
                    retcode=2; % nans in jacobian
                    return
                end
            end
        end
    end
end

