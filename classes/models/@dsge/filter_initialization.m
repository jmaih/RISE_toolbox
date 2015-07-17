function [init,retcode]=filter_initialization(obj,varargin)%T,R,steadystate,risk,transition_function,options
% filter_initialization - Initial conditions for filtering
%
% Syntax
% -------
% ::
%
%   [init,retcode]=filter_initialization(obj)
%
% Inputs
% -------
%
% - **T** [cell]: n x n elements representing the autoregressive part of the
%   first-order approximation of the model
%
% - **R** [cell]: n x ne elements representing the shock impact of the
%   first-order approximation of the model
%
% - **steadystate** [cell]: steady state of the model
%
% - **risk** [cell]: impact on the perturbation coefficient
%
% - **transition_function** [function handle]: function generating a
%   transition matrix for a given vector of endogenous
%
% - **options** [struct]: options with the most relevant beeing:
%
%   - **kf_ergodic** [{true}|false]: initialization at the ergodic
%       distribution
%
%   - **kf_init_variance** [{[]}|scalar]: initial variance factor (Harvey
%       scale factor). If not empty, the information in T and R is ignored.
%
%   - **kf_presample** [{0}|integer]: Number of observations to discard
%       before computing the likelihood.
%
%   - **kf_user_init** [{[]}|cell]: User-defined initialization. When not
%   empty, it can take three forms. {a0}, {a0,cov_a0}, {a0,cov_a0,PAI00}
%   where a0 is the initial state vector with the same order as the rows of
%   T, cov_a0 is the initial covariance of the state vector (same order as
%   a0) and PAI00 is the initial vector of regime probabilities.
%
%   - **kf_user_algo** [{''}|char|function handle]: User-defined filtering
%   algorithm. It should have the same inputs and outputs as e.g.
%   switching_divided_difference_filter.m. 
%
%   - **kf_householder_chol** [{false}|true]: if true, return the cholesky
%   decomposition when taking the householder transformation. This option
%   is primarily used in the switching divided difference filter.
%
% Outputs
% --------
%
% - **init** [struct]: initial condition of the Kalman filter
%
% - **retcode** [scalar]: 0 if there is no problem
%
% More About
% ------------
%
% The initial covariance is always computed using the current impact matrix
% of the shocks even when the anticipation horizon of the agents is greater
% than zero.
%
% Examples
% ---------
%
% See also: 


% diffuse initialization for all elements in the state vector including
% stationary ones. This is what Waggoner and Zha do, but then they take
% a presample. The intuition, I guess, is that the filter eventually
% updates everything to the correct values. In some other cases, one
% may want set the presample to the number of unit roots as I have seen
% some place before... the drawback is that if the model has lots of
% unit roots and the sample is short...

if isempty(obj)
    if nargout>1
        error([mfilename,':: with no input argument, the number of output arguments cannot exceed 1'])
    end
    defaults=struct(...
        'kf_ergodic',true,... 
        'kf_init_variance',[],... 
        'kf_presample',0,...
        'kf_user_init',[],...
        'kf_algorithm','lwz',...%     alternative is kn (Kim and Nelson)
        'kf_tol',1e-20,...
        'kf_filtering_level',3,...
        'kf_riccati_tol',1e-6,...
        'kf_nsteps',1,...
        'kf_user_algo','',...
        'kf_nan_future_obs_means_missing',true,...
        'kf_householder_chol',false);
    lyap_options=lyapunov_equation();
    defaults=utils.miscellaneous.mergestructures(defaults,lyap_options);
    init=defaults;
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
        warning('It is a good practive to have a presample period to initialize the filter under diffuse conditions')
    end
    init.start=kf_presample+1;
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
        horizon=max(obj.exogenous.shock_horizon)+1;
        h=obj.markov_chains.regimes_number;
        % covariance matrix of stochastic shocks
        %----------------------------------------
        SIGeta=repmat({eye(stoch_exo_nbr*horizon)},1,h);
        
        % load solution in order_var form
        %---------------------------------
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
        [grand_state]=Minimum_state_variable_filtration();
        
        Qfunc=prepare_transition_routine(obj);
        Qfunc=small_memo(Qfunc,iov,grand_state);
        % load the restrictions if any:
        %------------------------------
        sep_compl=complementarity_memoizer(obj);
        sep_compl=small_memo(sep_compl,iov,grand_state);
        
        % find the anticipated shocks
        anticipated_shocks=[];
        if ~isempty(sep_compl)
            anticipated_shocks=obj.exogenous.shock_horizon>0;
        end
        
        % consider only the first-order solution for the inialization of
        % state and covariance of state
        %-----------------------------------------------------------------
        Tx=T(1,:);
        nstates=numel(state_vars_location);
        Te=Tx;Te_det=Tx;Tsig=Tx;
        for rt=1:h
            Tsig{rt}=Tx{rt}(:,nstates+1);
            % separate deterministic from stochastic shocks
            %------------------------------------------------
            Te{rt}=reshape(Tx{rt}(:,nstates+2:end),[endo_nbr,exo_nbr,horizon]);
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
            a0=ss_star;
            if any(Tsig_star)
                Ix=eye(endo_nbr);
                Ix(:,state_vars_location)=Ix(:,state_vars_location)-Tx_star;
                a0=a0+Ix\Tsig_star;
            end
            if a0_given
                % correct the previous entries.
                if ~isstruct(kf_user_init{1})
                    error('the first entry of kf_user_init should be a structure')
                end
                endo_names=get(obj,'endo_list');
                % re-order and trim
                endo_names=endo_names(grand_state(iov));
                present_a0=false(1,endo_nbr);% this is the updated endo_nbr
                for iii=1:numel(endo_names)
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
                LxTx=Tx_star(state_vars_location,:);
                LxR=Te_star(state_vars_location,:);
                [vx,retcode]=lyapunov_equation(LxTx,LxR*LxR',obj.options);
                if ~retcode
                    P0=Tx_star*vx*Tx_star.'+Te_star*Te_star.';
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
                P0=repmat({P0},1,h);
                a0=repmat({a0},1,h);
            end
        end
        init=struct('a',{a0},'P',{P0},'PAI00',PAI00,'T',{T},...
            'Tx',{Tx},'Tsig',{Tsig},'Te',{Te},'Te_Te_prime',{Te_Te_prime},...
            'Te_det',{Te_det},'steady_state',{steady_state},...
            'new_order',new_order,'state_vars_location',state_vars_location,...
            'Qfunc',Qfunc,'anticipated_shocks',anticipated_shocks,...
            'sep_compl',sep_compl,'H',{obj.solution.H},'obs_id',obs_id,...
            'SIGeta',{SIGeta},'is_det_shock',is_det_shock,'horizon',horizon);
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
        
        function [grand_state]=Minimum_state_variable_filtration()
            grand_state=true(1,endo_nbr);
            if obj.options.kf_filtering_level==0
                grand_state=false(1,endo_nbr);
                affect_trans_probs=obj.endogenous.is_affect_trans_probs;
                grand_state(state_vars_location)=true;
                grand_state(obs_id)=true;
                grand_state(affect_trans_probs)=true;
                for ii=1:obj.options.solve_order
                    for reg=1:h
                        T{ii,reg}=T{ii,reg}(grand_state,:);
                        if ii==1
                            steady_state{reg}=steady_state{reg}(grand_state);
                        end
                    end
                end
                csgs=cumsum(grand_state);
                state_vars_location=csgs(state_vars_location);
                obs_id=csgs(obs_id);
                
                % update the number of endogenous
                %--------------------------------
                endo_nbr=sum(grand_state);
            end
        end
    end

end

function ffunc=small_memo(gfunc,iov,grand_state)
% sum functions are written to accept the endogenous variables in the order
% of obj.endogenous.name. This function makes it possible to have the
% variables in the inv_order_var order and still apply the functions of
% interest.
if isempty(gfunc)
    ffunc=gfunc;
else
    if nargin<2
        grand_state=[];
    end
    if isempty(grand_state)
        grand_state=1:numel(iov);
    end
    nx=numel(iov);
    ffunc=@do_it;
end

    function varargout=do_it(y)
        % the y vector could be longer if we are storing the filters, in
        % which case it is augmented with shocks at the end.
        x=zeros(nx,1);
        nxy=min(nx,numel(y));
        x(grand_state)=y(1:nxy,1);
        [varargout{1:nargout}]=gfunc(x(iov));
    end
end

