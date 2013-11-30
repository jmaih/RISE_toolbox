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
    obj=mergestructures(dsge_lc_solve(),msre_solve(),...
        struct('solver','msre_klein',...
        'solve_order',1,...
        'solve_expect_order',1,...
        'check_stability',true,...
        'vectorized_code_threshold',150));%'accelerate_solver',true
    return
end

if ~isempty(varargin)
    obj=set_options(obj,varargin{:});
    %     obj.options=mysetfield(obj.options,varargin{:});
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

structural_matrices=[];
% partitions=[];
structure_only=nargout>2;
solve_order=obj.options.solve_order;
params=obj.parameter_values;
if obj.is_optimal_policy_model || obj.estimation_under_way
    solve_order=1;
end
resolve_it=check_whether_to_resolve();
h=obj.markov_chains.regimes_number;

% solve zeroth order or steady state ... and measurement errors
%--------------------------------------------------------------
retcode=solve_zeroth_order();
if solve_order>0 && ~retcode
    % the parameters may have changed during the computation of the zeroth
    % order or steady state. That is why they are reloaded here
    params=obj.parameter_values;
%     eqtns_number=obj.equations.number;
    retcode=solve_first_order();
    if solve_order>1 && ~retcode
        retcode=solve_second_order();
        if solve_order>2
            error('orders greater than 2 not yet implemented')
        end
    end
end
if retcode && obj.options.debug
    decipher_error(retcode)
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
        solve_expect_order=max(1,max(obj.options.irf_horizon,obj.options.solve_expect_order));
        forward_expansion_has_changed=~isequal(obj.current_expansion_order,solve_expect_order);
        solution_algorithm_has_changed=~isequal(obj.current_solution_algorithm,obj.options.solver);
        solver_order_has_changed=~isequal(obj.current_solve_order,solve_order);
        resolve_it=obj.estimation_under_way ||...
            (~isempty(params) && ~isequal(obj.current_solution_parameters,params))||... % ~isequal(obj.current_solution_parameters,params)||...
            forward_expansion_has_changed||...
            solution_algorithm_has_changed||...
            strcmp(obj.options.solve_initialization,'random');
        
        if forward_expansion_has_changed
            obj.current_expansion_order=solve_expect_order;
        end
        if solution_algorithm_has_changed
            obj.current_solution_algorithm=obj.options.solver;
        end
        if solver_order_has_changed
            obj.current_solve_order=solve_order;
        end
    end

    function retcode=solve_zeroth_order()
        retcode=0;
        if resolve_it
            obj.solution=rise.initialize_solution_or_structure('solution',h);
            [obj,retcode]=compute_steady_state(obj);
            if ~retcode
                % evaluate steady-state transition matrix
                %----------------------------------------
                s0=1;s1=1;
                def=obj.solution.definitions{1};
                ss=obj.solution.ss{1};
                y=ss;
                x=zeros(sum(obj.exogenous.number),1);
                [obj.solution.Q,retcode]= online_function_evaluator(obj.func_handles.transition_matrix,y,x,ss,params(:,1),def,s0,s1);
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
    end
    function retcode=solve_first_order()
        retcode=0;
        resolve_it= resolve_it||isempty(obj.solution.m_x{1});
        if ~(resolve_it||structure_only)
            return
        end
        % evaluate derivatives up to the desired order
        %---------------------------------------------
        structural_matrices=rise.initialize_solution_or_structure('structure',h);
        retcode=evaluate_derivatives();
        if retcode
            return
        end
        % use the derivatives G1 to build Gp, Gc, Gm, Ge and Gt
        %------------------------------------------------------
        partitions=obj.model_derivatives.Endogenous_Shocks_Parameters(1).partitions;
        myfields=fieldnames(partitions);
        for s0=1:h
            for s1=1:h
                for ifield=1:numel(myfields)
                    ff=myfields{ifield};
                    structural_matrices.(['G',ff]){s0,s1}=sparse(obj.G1{s0,s1}(:,partitions.(ff)));
                end
            end
        end
        % add the planner information if necessary
        %-----------------------------------------
        
        if resolve_it
            % solve for m_x, m_x, m_sig
            %--------------------------
            [obj,retcode]=first_order_solver(obj,structural_matrices);
        end
    end
    function retcode=solve_second_order()
        retcode=0;
         resolve_it=resolve_it||isempty(obj.solution.m_x_x{1});
        if ~(resolve_it||structure_only)
            return
        end
        % use the derivatives G2 to build Gpp, Gpc, ..., Gtt
        %---------------------------------------------------
        partitions=obj.model_derivatives.Endogenous_Shocks_Parameters(2).partitions;
        myfields=fieldnames(partitions);
        for s0=1:h
            for s1=1:h
                for ifield=1:numel(myfields)
                    ff=myfields{ifield};
                    structural_matrices.(['G',ff]){s0,s1}=sparse(obj.G2{s0,s1}(:,partitions.(ff)));
                end
            end
        end
        if resolve_it
            % solve for m_x_x, m_x_e, m_x_sig, m_e_e, m_e_sig, m_sig_sig
            %-----------------------------------------------------------
            [obj,retcode]=second_order_solver(obj,structural_matrices);
        end
    end
    function retcode=evaluate_derivatives()
        %anonymous function for checking validity
        %----------------------------------------
        valid=@(x)~any(isnan(x(:))) && ~any(isinf(x(:))); % nans in jacobian
        % load the derivatives
        %---------------------
        derivative_type=obj.options.derivatives;
        switch derivative_type
            case 'symbolic'
                derivatives=obj.model_derivatives.Endogenous_Shocks_Parameters;%<-- derivatives=obj.model_derivatives.Endogenous_Shocks_Parameters{1};
            case {'numerical','automatic'}
                derivatives={obj.func_handles.vectorized_dynamic;
                    obj.func_handles.vectorized_dynamic_params};
        end
        xss=zeros(sum(obj.exogenous.number),1);
        Lead_lag_incidence=obj.Lead_lag_incidence;
        endo_nbr=obj.endogenous.number(2);
        
        the_leads=find(Lead_lag_incidence(:,1)>0);
        the_current=(1:endo_nbr)';
        the_lags=find(Lead_lag_incidence(:,3)>0);
        indices=[the_leads;the_current;the_lags];

        y=zeros(numel(indices),h);
        ss=zeros(endo_nbr,h);
        def=obj.solution.definitions;
        obj.G1=cell(h);
        obj.G2=cell(h);
        retcode=0;
        % first and second-order derivatives with respect to endogenous,
        % exogenous and switching parameters
        %---------------------------------------------------------------
        for s0=1:h
            for s1=1:h
                if s0==1
                    ss(:,s1)=obj.solution.ss{s1};
                    y(:,s1)=ss(indices,s1);
                    % y needs to be modified for the variables that have a
                    % nonzero balanced-growth path
                    if ~obj.is_stationary_model
                        if s1==1
                            [bgp_coefs,change_loc]=balanced_growth_path_powers();
                        end
                        y(change_loc,s1)=y(change_loc,s1).*obj.solution.bgp{s1}(indices(change_loc)).^bgp_coefs;
                    end
                end
                % Note: G(s0,s1) =: ps0(s0,s1)*F(s0)
                if solve_order>1
                    [obj.G1{s0,s1},obj.G2{s0,s1}]=online_function_evaluator(derivatives,y(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                    if ~valid(obj.G1{s0,s1})||~valid(obj.G2{s0,s1})
                        retcode=2; % nans in jacobian
                        return
                    end
                else
                    switch derivative_type
                        case 'symbolic'
                            [obj.G1{s0,s1}]=online_function_evaluator(derivatives,y(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                        case 'numerical'
                            [obj.G1{s0,s1}]=my_numerical_derivatives(derivatives,y(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                        case 'automatic'
                            [obj.G1{s0,s1}]=my_automatic_derivatives(derivatives,y(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                    end
                    if ~valid(obj.G1{s0,s1})
                        retcode=2; % nans in jacobian
                        return
                    end
                    % this is where I should put numerical and automatic
                    % derivatives as a choice. But then, this means that
                    % they should also allow the same input arguments as
                    % their symbolic counterpart.
                end
            end
        end
        % first and second-order derivatives of the planner objective
        %------------------------------------------------------------
        if obj.is_optimal_policy_model || obj.is_optimal_simple_rule_model
            [welfare_commitment_discount,first_derivatives,short_hessian]=...
                online_function_evaluator(obj.planner_system.LossComDiscHessJac,...
                ss(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s0); %#ok<ASGLU>
            structural_matrices.planner.objective=welfare_commitment_discount(1);
            structural_matrices.planner.commitment=welfare_commitment_discount(2);
            structural_matrices.planner.discount=welfare_commitment_discount(3);
            partitions=obj.planner_system.LossComDiscHessJac(3).partitions.cc;
            ncols=numel(partitions);
            weights=short_hessian(obj.planner_system.LossComDiscHessJac(3).partitions.cc);
            structural_matrices.planner.weights=reshape(weights,sqrt(ncols),sqrt(ncols));
        end
        function [junk,J]=my_numerical_derivatives(derivatives,y,xss,ss_i,param_i,def_i,s0,s1)
            junk=[];
            z=[y;xss];
            myfun=@(xx)online_function_evaluator(derivatives{1},xx,param_i,ss_i,def_i,s0,s1);
            J=rise_numjac(myfun,z);
            if ~isempty(switching)
                % now the z vector has to be a matrix in order to avoid a
                % breakdown occurring when no parameter is present in some
                % equations.
                z=z(:,ones(1,size(param_i,1)));
                myfun=@(xx)online_function_evaluator(derivatives{2},xx,z,ss_i,s0,s1);
                JP=rise_numjac(myfun,param_i);
                J=[J,JP];
            end
        end
        function  [junk,J]=my_automatic_derivatives(derivatives,y,xss,ss_i,param_i,def_i,s0,s1)
            junk=[];
            z=[y;xss];
            zz=rise_nad(z);
            J=online_function_evaluator(derivatives{1},zz,param_i,ss_i,def_i,s0,s1);
            if ~isempty(switching)
                JP=online_function_evaluator(derivatives{2},rise_nad(param_i),z,ss_i,s0,s1);
                J=[J,JP];
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

