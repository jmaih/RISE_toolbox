function [init,retcode]=filter_initialization(obj,varargin)
% FILTER_INITIALIZATION - Initial conditions for filtering
%
% ::
%
%
%   [init,retcode]=FILTER_INITIALIZATION(obj)
%
% Args:
%
%    - **obj** [rise|dsge]: model object
%
%    - **varargin** [name,value]: valid pairwise options with the most
%    relevant beeing:
%
%      - **kf_ergodic** [{true}|false]: initialization at the ergodic
%          distribution
%
%      - **kf_init_variance** [{[]}|scalar]: initial variance factor (Harvey
%          scale factor). If not empty, the information in T and R is ignored.
%
%      - **kf_presample** [{0}|integer]: Number of observations to discard
%          before computing the likelihood.
%
%      - **kf_filtering_level** [0|1|2|{3}]: 0: Likelihood only, 1: 0+filtered
%      series, 2: 1+ updated series, 3: 2+ smoothed series
%
%      - **kf_user_init** [{[]}|cell]: User-defined initialization. When not
%      empty, it can take three forms. {a0}, {a0,cov_a0}, {a0,cov_a0,PAI00}
%      where a0 is the initial state vector with the same order as the rows of
%      T, cov_a0 is the initial covariance of the state vector (same order as
%      a0) and PAI00 is the initial vector of regime probabilities.
%
%      - **kf_user_algo** [{''}|char|function handle]: User-defined filtering
%      algorithm. It should have the same inputs and outputs as e.g.
%      switching_divided_difference_filter.m.
%
%      - **kf_householder_chol** [{false}|true]: if true, return the cholesky
%      decomposition when taking the householder transformation. This option
%      is primarily used in the switching divided difference filter.
%
% Returns:
%    :
%
%    - **init** [struct]: initial conditions of the filter
%
%    - **retcode** [scalar]: 0 if there is no problem
%
% Note:
%
%    - The initial covariance is always computed using the current impact
%    matrix of the shocks even when the anticipation horizon of the agents is
%    greater than zero.
%
%    - use option "lyapunov_diffuse_factor" to set the variance for
%    nonstationary variables separately. If this is not the case the variance
%    for those variables will be 1 by default.
%
% Example:
%
%    See also: DSGE/FILTER


% diffuse initialization for all elements in the state vector including
% stationary ones. This is what Waggoner and Zha do, but then they take
% a presample. The intuition, I guess, is that the filter eventually
% updates everything to the correct values. In some other cases, one
% may want set the presample to the number of unit roots as I have seen
% some place before... the drawback is that if the model has lots of
% unit roots and the sample is short...

if isempty(obj)
    
    mydefaults=the_defaults();
    
    mydefaults=[mydefaults
        lyapunov_equation()];
    
    if nargout
        
        init=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if ~isempty(varargin)
    
    obj=set(obj,varargin{:});
    
end

kf_ergodic=obj.options.kf_ergodic;

kf_init_variance=obj.options.kf_init_variance;

kf_presample=obj.options.kf_presample;

kf_user_init=obj.options.kf_user_init;

kf_diffuse_all=~isempty(kf_init_variance);

[init,retcode]=initialize_filter();

if ~retcode
    
    kf_presample=max(kf_presample,0);
    
    if kf_diffuse_all && kf_presample==0
        
        warning('Diffuse conditions detected with zero presample')
        
        warning('It is a good practice to have a presample period to initialize the filter under diffuse conditions')
        
    end
    
    init.start=kf_presample+1;
    
    % re-inflate everything in order to store filters if necessary
    %-------------------------------------------------------------
    init=re_inflator(init,obj.options.kf_filtering_level);
    
end

%--------------------------------------------------------------------------
    function [init,retcode]=initialize_filter()
        
        if isempty(obj.solution)
            
            error('Model has not been solved')
            
        end
        
        stoch_exo_nbr=obj.exogenous.number(1);
        
        det_exo_nbr=obj.exogenous.number(2);
        
        exo_nbr=stoch_exo_nbr+det_exo_nbr;
        
        is_det_shock=obj.exogenous.is_observed;
        
        k=max(obj.exogenous.shock_horizon,[],2);
        
        horizon=max(k)+1;
        
        h=obj.markov_chains.regimes_number;
        
        % covariance matrix of stochastic shocks
        %----------------------------------------
        SIGeta=repmat({eye(stoch_exo_nbr*horizon)},1,h);
        
        % load solution in order_var form: the steady state is already
        % logged where necessary and the Qfunc is also in sync with that
        %------------------------------------------------------------------
        [T,~,steady_state,new_order,state_vars_location]=load_solution(obj,'ov');
        
        endo_nbr=obj.endogenous.number;
        
        iov(new_order)=1:numel(new_order);
        
        obs_id=real(obj.observables.state_id);
        
        % remove exogenous
        obs_id(obs_id==0)=[];
        
        % turn into the correct order
        obs_id=iov(obs_id);
        
        % If possible, trim solution for Minimum state filtering
        %---------------------------------------------------------
        [grand_order_var_state]=Minimum_state_variable_filtration();
        
        Qfunc=prepare_transition_routine(obj);
        
        Qfunc=small_memo(Qfunc,iov,grand_order_var_state);
        
        % load the restrictions if any:
        %------------------------------
        sep_compl=complementarity_memoizer(obj);
        
        sep_compl=small_memo(sep_compl,iov,grand_order_var_state);
        
        % find the anticipated shocks
        anticipated_shocks=[];
        
        if ~isempty(sep_compl)
            
            anticipated_shocks=any(obj.exogenous.shock_horizon>0,1);
            
        end
        
        % consider only the first-order solution for the inialization of
        % state and covariance of state
        %-----------------------------------------------------------------
        Tx=T(1,:);
        
        nstates=numel(state_vars_location);
        
        Te=Tx;Te_det=Tx;Tsig=Tx;
        
        for rt=1:h
            
            Tsig{rt}=Tx{rt}(:,nstates+1);
            % separate deterministic from stochastic shocks:
            % full is needed since the results are stored as sparse and we
            % are reshaping things in 3 dimensions
            % results are stored as sparse and so reshaping in
            % multidimensional requires full-ing
            %------------------------------------------------
            Te{rt}=reshape(full(Tx{rt}(:,nstates+2:end)),...
                [endo_nbr,exo_nbr,horizon]);
            
            Te_det{rt}=Te{rt}(:,is_det_shock,:);
            
            Te{rt}(:,is_det_shock,:)=[];
            
            % impact of endogenous state variables
            %--------------------------------------
            Tx{rt}=Tx{rt}(:,1:nstates);
            
        end
        
        a0_given=numel(kf_user_init)>0 && ~isempty(kf_user_init{1});
        
        P0_given=a0_given && numel(kf_user_init)>1 && ~isempty(kf_user_init{2});
        
        PAI00_given=numel(kf_user_init)>2 && ~isempty(kf_user_init{3});
        
        retcode=0;
        
        if PAI00_given
            
            PAI00=kf_user_init{3};
            
        else
            
            Q=Qfunc(steady_state{1});
            
            [PAI00,retcode]=initial_markov_distribution(Q,kf_ergodic);
            
        end
        
        Te_Te_prime=cell(h,1);
        
        [Tx_star,Tsig_star,Te_star,ss_star]=aggregate();
        
        P0=[];
        
        a0=[];
        
        if ~retcode
            
            stat_ind=stationary_index(obj);
            
            stat_ind=stat_ind(new_order);
            
            stat_ind=stat_ind(grand_order_var_state);
            
            is_state_var=false(size(stat_ind));
            
            is_state_var(state_vars_location)=true;
            
            is_stationary_state_var=stat_ind & is_state_var;
            
            % take the real part in order to avoid mixing up with growth
            %-------------------------------------------------------------
            Tsig_star=real(Tsig_star);
            
            a0=ss_star;
            
            if any(abs(Tsig_star)>1e-7)
                
                Ix=eye(endo_nbr);
                
                Ix(:,state_vars_location)=Ix(:,state_vars_location)-Tx_star;
                
                a0=a0+pinv(Ix)*Tsig_star;
                
            end
            
            if a0_given
                % correct the previous entries.
                if ~isstruct(kf_user_init{1})
                    
                    error('the first entry of kf_user_init should be a structure')
                    
                end
                
                endo_names=get(obj,'endo_list');
                
                % re-order and trim
                endo_names=endo_names(new_order);
                
                endo_names=endo_names(grand_order_var_state);
                
                present_a0=false(1,endo_nbr);% this is the updated endo_nbr
                
                for iii=1:numel(endo_names)
                    
                    name=endo_names{iii};
                    
                    if isfield(kf_user_init{1},name)
                        
                        if ~isempty(kf_user_init{1}.(name))
                            % we allow this to be empty if we just have the
                            % covariance information to come later...
                            a0(iii)=kf_user_init{1}.(name);
                            
                        end
                        
                        present_a0(iii)=true;
                        
                    end
                    
                end
                
            end
            
            if kf_diffuse_all
                
                P0=kf_init_variance*eye(endo_nbr);
                
            else
                
                P0=diag(obj.options.lyapunov_diffuse_factor*ones(endo_nbr,1));
                
                rows_ind=is_stationary_state_var;
                
                cols_ind=is_stationary_state_var(state_vars_location);
                
                LxTx=Tx_star(rows_ind,cols_ind);
                
                LxR=Te_star(rows_ind,:);
                
                [vx,retcode]=lyapunov_equation(LxTx,LxR*LxR',obj.options);
                
                if ~retcode
                    
                    LxTx=Tx_star(stat_ind,cols_ind);
                    
                    LxR=Te_star(stat_ind,:);
                    
                    P0(stat_ind,stat_ind)=LxTx*vx*LxTx.'+LxR*LxR.';
                    
                    P0=utils.cov.symmetrize(P0);
                    
                end
                
            end
            
            if P0_given
                
                np0=sum(present_a0);
                
                if issvector(kf_user_init{2})
                    
                    kf_user_init{2}=diag(kf_user_init{2});
                    
                end
                
                [nr,nc]=size(kf_user_init{2});
                
                if nr~=np0||nc~=np0
                    
                    error('second entry of kf_user_init should be consistent with first entry')
                    
                end
                
                P0(present_a0,present_a0)=kf_user_init{2};
                
            end
            
            if ~retcode
                
                if ~all(isfinite(a0))
                    
                    retcode=3002;
                    
                end
                
                if ~all(isfinite(P0(:)))
                    
                    retcode=30002;
                    
                end
                
                P0=repmat({P0},1,h);
                
                a0=repmat({a0},1,h);
                
            end
            
        end
        
        is_log_var=obj.log_vars;
        
        init=struct('a',{a0},'P',{P0},'PAI00',PAI00,'T',{T},...
            'Tx',{Tx},'Tsig',{Tsig},'Te',{Te},'Te_Te_prime',{Te_Te_prime},...
            'Te_det',{Te_det},'steady_state',{steady_state},...
            'new_order',new_order,'state_vars_location',state_vars_location,...
            'Qfunc',Qfunc,'anticipated_shocks',anticipated_shocks,...
            'sep_compl',sep_compl,'H',{obj.solution.H},'obs_id',obs_id,...
            'SIGeta',{SIGeta},'is_det_shock',is_det_shock,'horizon',horizon,...
            'k',k,'is_log_var_new_order',is_log_var(new_order));
        
        function [Tx_star,Tsig_star,Te_star,ss_star]=aggregate()
            
            Tx_star=0;
            
            Te_star=0;
            
            Tsig_star=zeros(endo_nbr,1);
            
            ss_star=0;
            
            for ireg=1:h
                
                Tx_star=Tx_star+PAI00(ireg)*Tx{ireg};
                
                Te_star=Te_star+PAI00(ireg)*Te{ireg}(:,:,1);
                
                Te_Te_prime{ireg}=Te{ireg}(:,:,1)*Te{ireg}(:,:,1)';
                
                Tsig_star=Tsig_star+PAI00(ireg)*Tsig{ireg};
                
                ss_star=ss_star+PAI00(ireg)*steady_state{ireg};
                
            end
            
        end
        
        function [grand_order_var_state]=Minimum_state_variable_filtration()
            
            grand_order_var_state=true(1,endo_nbr);
            
            if obj.options.kf_filtering_level==0
                
                grand_order_var_state=false(1,endo_nbr);
                
                grand_order_var_state(state_vars_location)=true;
                
                grand_order_var_state(obs_id)=true;
                
                % the variables affecting the transition probabilities need
                % to be re-ordered before locating them in grand_state
                affect_trans_probs=obj.endogenous.is_affect_trans_probs;
                
                grand_order_var_state(affect_trans_probs(new_order))=true;
                
                for ii=1:obj.options.solve_order
                    
                    for reg=1:h
                        
                        T{ii,reg}=T{ii,reg}(grand_order_var_state,:);
                        
                        if ii==1
                            
                            steady_state{reg}=steady_state{reg}(grand_order_var_state);
                            
                        end
                        
                    end
                    
                end
                
                csgs=cumsum(grand_order_var_state);
                
                state_vars_location=csgs(state_vars_location);
                
                obs_id=csgs(obs_id);
                
                % update the number of endogenous
                %--------------------------------
                endo_nbr=sum(grand_order_var_state);
                
            end
            
        end
        
    end

end

function ffunc=small_memo(gfunc,iov,grand_order_var_state)
% some functions are written to accept the endogenous variables in the
% order of obj.endogenous.name. This function makes it possible to have the
% variables in the inv_order_var (alphabetical) order and still apply the
% functions of interest.
if isempty(gfunc)
    
    ffunc=gfunc;
    
else
    
    if nargin<2
        
        grand_order_var_state=[];
        
    end
    
    if isempty(grand_order_var_state)
        
        grand_order_var_state=1:numel(iov);
        
    end
    
    nx=numel(iov);
    
    ffunc=@do_it;
    
end

    function varargout=do_it(y)
        % the y vector could be longer if we are storing the filters, in
        % which case it is augmented with shocks at the end.
        x=zeros(nx,1);
        
        nxy=min(nx,numel(y));
        
        x(grand_order_var_state)=y(1:nxy,1);
        
        [varargout{1:nargout}]=gfunc(x(iov));
        
    end

end

function mydefaults=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

mydefaults={
    'kf_ergodic',true,@(x)islogical(x),'kf_ergodic must be true or false'
    
    'kf_init_variance',[],...
    @(x)isempty(x)||(num_fin(x) && x>0),...
    'kf_init_variance must be a positive scalar'
    
    'kf_presample',0,@(x)num_fin_int(x),...
    'kf_presample must be zero or a positive integer'
    
    'kf_user_init',[],@(x)isempty(x)||(iscell(x) && numel(x)<=3),...
    'kf_user_init must be a cell with number of entries <=3'
    
    'kf_tol',1e-20,@(x)num_fin(x) && x>=0,'kf_tol must be >=0'
    
    'kf_filtering_level',3,@(x)any(x==(0:3)),...
    'kf_filtering_level must be 0,1,2 or 3'
    
    'kf_riccati_tol',1e-6,@(x)num_fin(x) && x>=0,...
    'kf_riccati_tol must be >=0'
    
    'kf_nsteps',1,@(x)num_fin_int(x) && x>=1,...
    'kf_nsteps must be an integer >=1'
    
    'kf_user_algo','',...
    @(x)isempty(x)||iscell(x)||ischar(x)||isa(x,'function_handle'),...
    'kf_user_algo must be a char, a cell, or a function handle'
    
    'kf_nan_future_obs_means_missing',true,...
    @(x)islogical(x),...
    'kf_nan_future_obs_means_missing must be a logical'
    
    'kf_householder_chol',false,...
    @(x)islogical(x),...
    'kf_householder_chol must be a logical'
    };

end

function init=re_inflator(init,kf_filtering_level)
m_orig=size(init.a{1},1);
if kf_filtering_level
    is_det_shock=init.is_det_shock;
    exo_nbr=numel(is_det_shock);
    horizon=init.horizon;
    nshocks=exo_nbr*horizon;
    h=numel(init.PAI00);
    nstates=numel(init.state_vars_location);
    for st=1:h
        init.a{st}=[init.a{st};zeros(nshocks,1)];
        init.steady_state{st}=[init.steady_state{st};zeros(nshocks,1)];
        init.P{st}=[init.P{st},zeros(m_orig,nshocks)
            zeros(nshocks,m_orig),eye(nshocks)];
        for io=1:size(init.T,1)
            ncols=size(init.T{io,st},2);
            batch=zeros(nshocks,ncols);
            if io==1
                % place the identities
                batch(:,nstates+2:end)=eye(nshocks);
            end
            init.T{io,st}=[init.T{io,st};batch];
            if io==1
                % collect the Te and Te_det instead of recomputing them
                %-------------------------------------------------------
                Te_all=reshape(full(init.T{io,st}(:,nstates+2:end)),...
                    [m_orig+nshocks,exo_nbr,horizon]);
                init.Te{st}=Te_all(:,~is_det_shock,:);
                init.Te_det{st}=Te_all(:,is_det_shock,:);
                
            end
        end
        % decoupled first order
        %----------------------
        init.Tx{st}=[init.Tx{st};zeros(nshocks,nstates)];
        init.Tsig{st}=[init.Tsig{st};zeros(nshocks,1)];
        
        Te_Te=init.Te{st}(:,:)*init.Te{st}(:,:).';
        Te_Te(1:m_orig,1:m_orig)=init.Te_Te_prime{st};
        init.Te_Te_prime{st}=Te_Te;
    end
end
init.m_orig=m_orig;
init.m=size(init.a{1},1);

end