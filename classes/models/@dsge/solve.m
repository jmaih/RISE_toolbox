function [obj,retcode,structural_matrices]=solve(obj,varargin)
% solve - solves dsge model
%
% ::
%
%
%   [obj,retcode,structural_matrices]=solve(obj)
%
%   [obj,retcode,structural_matrices]=solve(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge]: scalar or vector of model objects. The optional
%    options below come in pairs.
%
%    - **solve_accelerate** [{false}|true]: Accelerate or do not accelerate
%    the solving
%
%    - **solve_check_stability** [{true}|false]: check stability of Markov
%    switching models while solving. The stability of constant-parameter
%    models is always checked whether that of markov-switching models is
%    optional. This is because (1) the procedure is computationally intensive
%    and (2) there is no define stability criterion under endogenous switching
%
%    - **solve_derivatives_type** [numeric|automatic|{symbolic}]: choice of
%    derivatives
%
%    - **solve_order** [integer|{1}]: order of approximation
%
%    - **solve_shock_horizon** [integer|{0}|struct|cell]: anticipation horizon
%    of shocks beyond the current period. When the input is :
%      - an integer : all the shocks receive the same anticipation horizon
%      - a struct : the fields are the names of the shocks whose horizon is to
%      be modified. e.g. struct('ea',4,'eb',3) means shock ea has horizon 4
%      while shock eb has horizon 3
%      - a cell : the cell must have two colums. Each row in the first column
%      holds the name of a particular shock and each row in the second column
%      holds the horizon of the shock. e.g. {'ea',4;'eb',3}
%
%    - **solve_alternatives_file2save2** [{[]}|char]: name of the file to
%    write the results of the alternative solutions
%
%    - **solve_alternatives_nsim** [integer|{100}]: Number of initial guesses
%    for the solution sampled randomly, when attempting to find all possible
%    solutions.
%
%     - **solve_function_mode** [{explicit/amateur}|vectorized/professional|disc]
%      - in the **amateur** or **explicit** mode the functions are kept in
%      cell arrays of anonymous functions and evaluated using for loops
%      - in the **vectorized** or **professional** mode the functions are
%      compacted into one long and unreadable function.
%      - in the **disc** mode the functions are written to disc in a
%      subdirectory called routines.
%
%    - **solve_initialization** [{'backward'}|'zeros'|'random']: Type of
%    initialization of the solution of switching models:
%      - **backward** : the initial guess is the solution of the model without
%      forward-looking terms
%      - **zeros** : the initial guess is zero
%      - **random** : the initial guess is random
%
%    - **solve_linsyst_user_algo** [{''}|cell|char]: user-defined solver for
%    linear systems. It should be possible to call the function as
%    [X,FLAG,RELRES,ITER,RESVEC] = bicgstabl(A,B,TOL,MAXIT,varargin). Hence,
%    the function should accept at least 4 inputs where:
%      - **A** : matrix or function handle
%      - **B** : right-hand side of system to solve
%      - **TOL** : tolerance level
%      - **MAXIT** : maximum number of iterations
%    The function must return 5 outputs as in Matlab's tfqmr, bicgstab,
%    bicgstabl, gmres, etc. :
%      - **X** : the solution
%      - **FLAG** : convergence flag
%      - **RELRES** : relative residual NORM(B-A*X)/NORM(B)
%      - **ITER** : iteration number
%      - **RESVEC** : residual vector for various iterations
%    As far as the outputs are concerned, only the two first are relevant for
%    RISE and so the three last can be empty. But they must be returned by the
%    function all the same.
%
%    - **solver** [{[]}|char|user_defined]: solver for the dsge model. The
%    following are possible:
%      - **loose_commitment** : RISE automatically uses this when it detects
%      an optimal policy model to be solved under commitment, discretion or
%      loose commitment
%      - **rise_1** : default solver for constant-parameter dsge models.
%      Similar to the dynare solver.
%      - **mfi** : functional interation solver: default for switching dsge
%      models
%      - **mnk** : newton solver with kronecker products
%      - **mn** : newton solver without kronecker products
%      - **mfi_full** : full version of mfi, that does not exploit sparsity
%      - **mnk_full** : full version of mnk, that does not exploit sparsity
%      - **mn_full** : full version of mn, that does not exploit sparsity
%      - **fwz** : Farmer, Waggoner and Zha (2011) solver
%      - **user_defined** : In this case the function must take as inputs:
%          - **Gplus01** [n x n x h x h array]: where the 3rd dimension
%          represent the current regime and the 4th dimension is the future
%          regime.
%          - **A0** [n x n x h array]: matrix of contemporaneous variables
%          - **Aminus** [n x n x h array]: matrix of backward-looking terms
%          - **Q** [square matrix]: transition matrix, in which the rows are
%          the current period and the columns the next period
%          - **T0** [n x h x h array]: initial guess of the solution
%          - **tol** [scalar]: A tolerance criterion
%          - **maxiter** [scalar]: the maximum number of iterations
%        The function should return as outputs
%          - **Tz_pb** [n x h x h array]: solution of the problem given the
%          inputs
%          - **eigenvalues** [empty|vector]: optional vector of eigenvalues of
%          the problem
%          - **retcode** : 0 if no problem found when solving the problem
%       The function can be passed along as a string, a function handle or a
%       cell array in which case the first element of the cell contains the
%       function itself and the remaining cells contain further arguments
%       required by the function but unknown to RISE.
%
%    - **solve_log_approx_vars** [char|cellstr|{[]}]: List of variables for
%    which we want to take a log expansion (x_t-x_ss)/x_ss, which approximates
%    log(x_t/x_ss) for x_t/x_ss close to 1.
%
%    - **solve_automatic_differentiator** [function_handle|@aplanar.diff|{@aplanar_.idff}]:
%    automatic differentiator engine
%
%    - **solve_higher_order_solver** ['dsge_solver_ha'|{'dsge_solver_h'}]: The
%    two solvers only differ in the way they handle memory. The default
%    solver (dsge_solver_h) keeps lots of data in memory while the
%    alternative one (dsge_solver_ha) recomputes some of the variables in
%    order to economize on memory. It is not clear which one performs the best
%
%    - **solve_occbin** [integer|{empty}]: Solves the occasionally
%    binding constraints model ala Kulish (2013) or Guerrieri and
%    Iacoviello(2015). It is then assumed that the transition matrix is
%    diagonal. if not empty, the value entered is the reference regime.
%
%    - **solve_bgp** [true|{false}]: Solves the model as non-stationary.
%
%    - **solve_sstate_blocks** [true|{false}]: blockwise solution of the
%    steady state.
%
%    - **solve_bgp_shift** [numeric|{5}]: shift/lead of the static model for
%    solving the balanced growth path.
%
%    - **solve_reuse_solution** [false|{true}]: re-use previous solution as
%    initial conditions
%
%    - **solve_linear** [true|{false}]: set this to true in order to use more
%    efficient algorithms for solving the steady state when the original model
%    is truly linear.
%
%    - **solve_time_consistent_switch** [true|{false}]: if true, at the
%    evaluation of the derivatives, leads are evaluated at the steady state of
%    the future regime
%
%    - **solve_derivatives_only** [true|{false}]: if true, the derivatives are
%    computed but the model is not solved. The derivatives can be collected in
%    the third output argument of the function.
%
%    - **steady_state_file** [char|function_handle|{''}]:
%
%    - **steady_state_use_steady_state_model** [false|{true}]:
%
%    - **steady_state_solver** [char|function_handle|{'lsqnonlin'}]:
%
%    - **steady_state_algorithm** [char|function_handle|{{'levenberg-marquardt',2*.005}}]:
%
%    - **steady_state_unique** [true|{false}]:
%
%    - **steady_state_imposed** [true|{false}]:
%
%    - **steady_state_loop** [true|{false}]:
%
%    - **steady_state_fixed** [true|{false}]:
%
%    - **steady_state_use_jacobian** [true|{false}]: In linear models, the
%    jacobian is always used. In the nonlinear case, however, the true
%    jacobian tends not to work as well as its finite differences
%    approximation. Hence by default, we do not invoke the use of the true
%    jacobian.
%
%    - **steady_state_endo_param_swap** [cell array|{}]: When not empty, it is
%    a cell array with n rows and 3 columns, where n is the number of
%    restrictions, the first column gathers the names of the endogenous
%    variables whose values are given in the second column, the third column
%    includes the parameters names that are endogenized.
%
% Returns:
%    :
%
%    - **obj** [rise|dsge]: scalar or vector of model objects
%
%    - **retcode** [integer]: if 0, no problem was found when solving the
%    model. If positive, the reason for the problem can be obtained by running
%    decipher(retcode)
%
%    - **structural_matrices** [struct]: Structure holding various types of
%    derivatives that are used in the computation of the solution
%
% Note:
%
%    - constant-parameter models can also be solved using the procedures for
%    switching parameter models. There is no guarantee, however, any obtained
%    solution will be determinate (unique and stable)
%
%    - determinacy is only checked for constant-parameter models. For
%    determinacy computations for switching dsge models with constant
%    transition matrices, see Cho(2014).
%
%    - In RISE we consider that it is in the nature of nonlinear models to
%    have many solutions. Obtaining even one of those solution is usually
%    difficult and a concept such as determinacy is NEVER used in the solving
%    of models with, say, global methods.
%
% Example:
%
%    See also:

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

if isempty(obj)
    
    mydefaults=the_defaults();
    
    mydefaults=[mydefaults
        fix_point_iterator()% add the defaults from fix point iterator
        optimal_policy_solver_h(obj)]; % loose commitment 
    
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj % to avoid displaying the empty object
        
        disp_defaults(mydefaults);
        
    end
    
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

% this flag forces resolve to get the structural matrices of the system
%-----------------------------------------------------------------------
is_struct_mat_out=nargout>2;

% initialize elements that will be used by different sub-functions
%-----------------------------------------------------------------
structural_matrices=[];

if ~isempty(obj.options.solve_occbin)
    
    obj.options.solve_order=min(obj.options.solve_order,1);
    
    obj.options.irf_type='irf';
    
end

solve_order=obj.options.solve_order;

resolve_it=check_whether_to_resolve();

h=obj.markov_chains.small_markov_chain_info.regimes_number;

% solve zeroth order or steady state ... and measurement errors
%--------------------------------------------------------------
retcode=solve_zeroth_order();

if ~retcode
    
    [structural_matrices,retcode]=dsge_tools.evaluate_all_derivatives(obj,structural_matrices);
    
    if obj.options.solve_derivatives_only
        
        return
        
    end
    
end

if solve_order>0 && ~retcode && resolve_it
    
    solve_occbin=obj.options.solve_occbin;
    
    if iscell(solve_occbin),solve_occbin=solve_occbin{1}; end
    
    obj.options.occbin.do_it=do_occbin(solve_occbin,...
        obj.markov_chains.regimes_number);
    
    if ~retcode
        
        loose_com_col=[];
        
        is_loose_commit=false;
        
        if obj.is_optimal_policy_model
            
            is_loose_commit=(structural_matrices.planner.commitment{1}<1||...
                strcmp(obj.options.solver,'loose_commitment'));
            
            loose_com_col=find(strcmp(obj.markov_chains.chain_names,parser.loose_commit));
            
            if (~isempty(loose_com_col)||is_loose_commit) && solve_order>1
                
                error('Loose commitment model only solved up to first order so far')
                
            end
            
        end
        
        if is_loose_commit
            % for loose commitment, only order 1 is available
            %-------------------------------------------------
            obj.options.solve_order=1;
            
            [T,eigval,retcode,obj]=optimal_policy_solver_h(obj,structural_matrices);
            
            % set the name of the solver
            %---------------------------
            if isempty(obj.options.solver)
                
                obj.options.solver='loose_commitment';
                
            end
            
        else
            
            if ischar(obj.options.solve_higher_order_solver)
                
                obj.options.solve_higher_order_solver=...
                    str2func(obj.options.solve_higher_order_solver);
                
            end
            
            [T,eigval,retcode,obj]=...
                obj.options.solve_higher_order_solver(obj,structural_matrices);
            
        end
        
        if ~retcode
            
            if obj.is_optimal_policy_model||obj.is_optimal_simple_rule_model
                
                planner=structural_matrices.planner;
                
                % change the ordering of the weight for the user. OSR will use
                % the structural matrices, which are ordered according to
                % order_var
                iov=obj.inv_order_var;
                
                for s00=1:h
                    
                    planner.weights{s00}=planner.weights{s00}(iov,iov);
                    
                end
                
                obj=set_planner_derivatives(obj,planner);
                
                clear planner
                
            end

            % expand solution to account for loose commitment
            %------------------------------------------------
            [obj,T]=expand_optimal_policy_solution(obj,T,loose_com_col);
            
            % Take the log approximation if required
            solve_log_approx_vars=obj.options.solve_log_approx_vars;
            
            if ~isempty(solve_log_approx_vars)
                
                retcode=log_expansion();
                
            end
            
            if ~retcode
                
                inv_order_var=obj.inv_order_var;
                
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
                
                if obj.options.solve_check_stability && ...
                        ~any(strcmp(obj.options.solver,{'rise_1','klein','aim','sims'}))%
                    
                    if ~obj.is_stable_system
                        
                        retcode=25; % system unstable
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

if retcode
    
    if obj.options.debug
        
        utils.error.decipher(retcode)
        
    end
    
else
    
    obj.warrant_resolving=false;
    
    obj.warrant_setup_change=false;
    
    % do the bvar-dsge setup
    do_bvar_dsge=obj.is_dsge_var_model && obj.options.dsgevar_var_regime;
    
    if do_bvar_dsge
        
        [obj,retcode]=bvar_dsge(obj);
        
    end
    
end

    function resolve_it=check_whether_to_resolve()
        % forward expansion order
        % if order==0, it means there is no ancitipation. It does not mean
        % that there no shocks. If order==0, this will later on help set
        % the conditional information to nan. But for the purpose of
        % expanding the model, the expansion order is at least 1.
        resolve_it=obj.estimation_under_way||obj.warrant_resolving||...
            strcmp(obj.options.solve_initialization,'random')||...
            ~isfield(obj.solution,'Tz')||is_struct_mat_out;
        
    end

    function retcode=solve_zeroth_order()
        
        retcode=0;
        
        if resolve_it
            % This may take a lot of space in large models or higher-order
            % approximations
            %--------------------------------------------------------------
            if ~isempty(obj.solution) && ...
                    isfield(obj.solution,'Tz') && ...
                    obj.options.solve_reuse_solution
                
                main_fields={'ss','bgp','Tz'};
                
                allfields=fieldnames(obj.solution);
                
                bad_fields=allfields-main_fields;
                
                obj.old_solution=rmfield(obj.solution,bad_fields);
                
                ov=obj.order_var;
                
                for ireg_=1:h
                    
                    obj.old_solution.Tz{ireg_}=obj.old_solution.Tz{ireg_}(ov,:);
                    
                end
                
            else
                % make sure the old solution is truly empty
                obj.old_solution=[];
                
            end
            
            if isempty(obj.solution)||~obj.estimation_under_way
                
                obj.solution=struct();%dsge.initialize_solution_or_structure('solution',h);
                
                obj.solution.H=cell(1,h);
                
            end
            
            % check whether I am not changing some data types in here to
            % explain why obj on the left-hand side is updated so slowly. I
            % could try different things :
            % 1- only output the solution and the updated parameters...
            % 2- put most of those things are sub-functions or nested
            % functions.
            [obj,structural_matrices,retcode]=compute_steady_state(obj);
            
            % put a copy of the old solution into the structural matrices
            %-------------------------------------------------------------
            structural_matrices.old_solution=obj.old_solution;
            
            % measurement errors
            %-------------------
            if ~retcode
                
                measure_flag=~isempty(obj.measurement_errors_restrictions);
                
                if measure_flag
                    
                    params=obj.parameter_values;
                    
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

    function retcode=log_expansion()
        
        retcode=0;
        
        if ischar(solve_log_approx_vars)
            
            solve_log_approx_vars=cellstr(solve_log_approx_vars);
            
        end
        
        log_vars=false(1,obj.endogenous.number);
        
        log_vars(locate_variables(solve_log_approx_vars,obj.endogenous.name))=true;
        
        % cannot be some original log_var
        %---------------------------------
        double_log=log_vars & obj.endogenous.is_log_var;
        
        if any(double_log)
            
            disp(obj.endogenous.name(double_log))
            
            error(['The variables above were already in log-form. ',...
                'A further log-expansion cannot be taken on them'])
            
        end
        
        % cannot be observable
        %----------------------
        if obj.observables.number(1)
            
            obs_id=real(obj.observables.state_id);
            
            log_obs=log_vars(obs_id);
            
            if any(log_obs)
                
                disp(obj.endogenous.name(obs_id(log_obs)))
                
                error('The variables above are observed. They cannot be log-expanded')
                
            end
            
        end
        
        tt=obj.locations.after_solve.t;
        
        pb_pos=[tt.p,tt.b];
        
        npb=numel(pb_pos);
        
        order_var=obj.order_var;
        
        nstates=size(T.Tz{1},2);
        
        ss_state_list=ones(1,nstates);
        
        for ireg_=1:h
            
            ss__=obj.solution.ss{ireg_}(:).';
            
            obj.solution.ss{ireg_}(log_vars)=0;
            
            % cannot have a zero steady state
            %---------------------------------
            if any(abs(ss__(log_vars))<obj.options.fix_point_TolFun)
                
                if obj.options.debug
                    
                    zero_sstate=abs(ss__)<obj.options.fix_point_TolFun & log_vars;
                    
                    disp(obj.endogenous.name(zero_sstate))
                    
                end
                
                retcode=27;
                
            end
            
            if retcode
                
                return
                
            end
            
            % Set to 1, the temporary steady state of the variables that
            % are not log expanded. Otherwise we will be forcing all the
            % variables to be log-expanded including those that have
            % steady state zero.
            ss__(~log_vars)=1;
            
            % now operate in the order_var mode
            %----------------------------------
            ss__=ss__(order_var);
            
            ss_state_list(1:npb)=ss__(pb_pos);
            
            zkz=1;
            
            Tz_='T';
            
            for io=1:obj.options.solve_order
                
                zkz=kron(zkz,ss_state_list);
                
                Tz_=[Tz_,'z']; %#ok<AGROW>
                
                T.(Tz_){ireg_}=bsxfun(@times,T.(Tz_){ireg_},zkz);
                
                % we need to full this because the various multiplications
                % transform to sparse and create problems with reshape
                % later on.
                T.(Tz_){ireg_}=full(bsxfun(@rdivide,T.(Tz_){ireg_},ss__(:)));
                
            end
            
        end
        
    end

end

function do_it=do_occbin(occbin,nregs)

do_it=@(x)true;

if ~isempty(occbin)
    
    if ~ismember(occbin,1:nregs)
        
        error('solve_occbin cannot exceed the number of regimes')
    
    end
    
    do_it=@(x)x==occbin;
    
end

end

function mydefaults=the_defaults()
% s: implies setup change
% r: implies re-solving of the model

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

functype=@(x)isempty(x)||ischar(x)||isa(x,'function_handle');

algotype=@(x)functype(x)||iscell(x);

    mydefaults={
        'solve_function_mode(sr)','explicit',...
        @(x)ismember(x,{'explicit','disc','vectorized'}),...
        'solve_function_mode must be one of "explicit","disc","vectorized"' 
        
        'solve_derivatives_type(sr)','symbolic', ...
        @(x)ismember(x,{'symbolic','numerical','automatic'}),...
        'solve_derivatives_type must be one of "symbolic","numerical","automatic"'
        
        'solve_bgp(sr)',false,@(x)islogical(x),'solve_bgp must be true or false'
        
        'solve_sstate_blocks(sr)',false,@(x)islogical(x),'solve_sstate_blocks must be true or false'
        
        'solve_linear(sr)',false,@(x)islogical(x),'solve_linear must be true or false'
        
        'solve_time_consistent_switch(r)',false,@(x)islogical(x),'solve_time_consistent_switch must be true or false'
        
        'steady_state_file(sr)','',@(x)functype(x),...
        'steady_state_file must be empty or char or a function handle'
        
        'steady_state_use_steady_state_model(sr)',true,@(x)islogical(x),...
        'steady_state_use_steady_state_model must be true or false'
        
        'steady_state_solver(sr)','lsqnonlin',...
        @(x)ischar(x)||isa(x,'function_handle'),...
        'steady_state_solver must be char or a function handle'
        
        'steady_state_algorithm(sr)',{'levenberg-marquardt',2*.005},...
        @(x)ischar(x)||iscell(x),'steady_state_algorithm must be a char or a cell'
        
        'steady_state_unique(sr)',false,@(x)islogical(x),...
        'steady_state_unique must be true or false'
        
        'steady_state_imposed(sr)',false,@(x)islogical(x),...
        'steady_state_imposed must be true or false'
        
        'steady_state_loop(sr)',false,@(x)islogical(x),...
        'steady_state_loop must be true or false'
        
        'steady_state_use_jacobian(sr)',false,@(x)islogical(x),...
        'steady_state_use_jacobian must be true or false'
        
        'steady_state_endo_param_swap(sr)',[],...
        @(x)isempty(x)||(iscell(x) && size(x,2)==3),...
        ['steady_state_endo_param_swap must be empty or a ',...
        'n x 3 cell {endogNames,endogValues,paramNames}']
        
        'steady_state_fixed(sr)',false,@(x)islogical(x),...
        'steady_state_fixed must be true or false'
        
        'solver(r)',[],@(x)algotype(x),...
        'solver can be empty, a char, a cell or a function handle'
        
        'solve_order(r)',1,@(x)num_fin_int(x) && x>0 && x<=5,...
        'solve_order must be a finite integer in [1,5]'
        
        'solve_shock_horizon(r)',[],@(x)isstruct(x)||num_fin_int(x),...
        'solve_shock_horizon must be >=0 or a structure whose fieldnames are the shocks'
        
        'solve_initialization(r)','backward',...
        @(x)ismember(x,{'zeros','backward','random'}),...
        'solve_initialization must be ''zeros'',''backward'' or ''random'''
        
        'solve_log_approx_vars(r)',[],...
        @(x)isempty(x)||ischar(x)||iscellstr(x),...
        'solve_log_approx_vars must be a char or a cellstr'
        
        'solve_accelerate(r)',false,@(x)islogical(x),...
        'solve_accelerate must be true or false'
        
        'solve_linsyst_user_algo(r)','',@(x)algotype(x),...
        'solve_linsyst_user_algo must be a char, a cell or a function handle'
        
        'solve_check_stability',true,@(x)islogical(x),...
        'solve_check_stability must be true or false'
        
        'solve_automatic_differentiator',@aplanar.diff,...
        @(x)isa(x,'function_handle'),...
        'solve_automatic_differentiator must be a function handle'
        
        'solve_higher_order_solver','dsge_solver_h',@(x)algotype(x),...
        'solve_higher_order_solver must be a char, a cell or a function handle'
        
        'solve_occbin',[],@(x)isempty(x)||...
        (iscell(x) && numel(x)==3 && (num_fin_int(x{1}) && x{1}>0))||...
        (num_fin_int(x) && x>0),...
        ['solve_occbin must be a positive integer or a three-element cell ',...
        'whose first element is a positive integer, the second a structure mapping ',...
        'the restrictions and the markov chains and the third ',...
        'the regimes']
        
        'solve_bgp_shift',5,@(x)(num_fin_int(x) && x>0),...
        'solve_bgp_shift must be a positive integer'
        
        'solve_reuse_solution',false,@(x)islogical(x),...
        'solve_reuse_solution must be true or false'
        
        'solve_derivatives_only',false,@(x)islogical(x),...
        'solve_derivatives_only must be true or false'
        };
        
end