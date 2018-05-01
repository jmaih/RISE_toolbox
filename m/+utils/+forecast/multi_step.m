function [sims,regimes,retcode,Qt,myshocks]=multi_step(y0,ss,T,state_vars_location,options)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% Qt below may be used to get the time series of the transition matrix
endo_nbr=size(y0.y,1);

PAI=options.PAI;

Qfunc=options.Qfunc;

% use first page only as it is the most relevant
regimes=vec(y0.rcond.data(:,:,1));

switch_rule=y0.switch_rule;

inv_order_var=y0.inv_order_var;

shocks=y0.econd.data;

cond_shocks_id=[];

if options.k_future
    
    if isfield(options,'shock_structure')
        
        cond_shocks_id=any(options.shock_structure,2);
        
    else
        
        cond_shocks_id=any(isnan(shocks(:,2:end,1)),2);
        
    end
    
end

forecast_conditional_hypothesis=options.forecast_conditional_hypothesis;

options=rmfield(options,{'PAI','Qfunc','y'});

[solve_order,h]=size(T);

Qt=[];

retcode=0;

penalty=1e+6;

span=options.nsteps+options.burn;

sims=nan(endo_nbr,options.nsteps);

nsv=numel(state_vars_location);

nx=(size(T{1,1},2)-nsv-1)/(options.k_future+1); %<---size(shocks,1);

condforkst=~isempty(y0.ycond.data)||...
    (...
    ~isempty(shocks) && (any(isnan(shocks(:)))||any(isinf(shocks(:)))||...
    max(max(abs(shocks(:,:,1)-shocks(:,:,2))))>1e-19||...
    max(max(abs(shocks(:,:,1)-shocks(:,:,3))))>1e-19||...
    max(max(abs(shocks(:,:,2)-shocks(:,:,3))))>1e-19)...
    )||...
    options.forecast_conditional_sampling_ndraws>1||...
    options.simul_frwrd_back_shoot;

if options.simul_anticipate_zero
    
    NotAnt_T=no_anticipation_solution();
    
end
%-------------------------------------
all_known_states=~any(isnan(regimes));

unique_state=numel(PAI)==1||(all_known_states && all(regimes==regimes(1)));
% model_is_linear=size(T,1)==1;
if unique_state
    
    if isnan(regimes(1))
        % no states were provided. Uniqueness in this case implies that
        % we have only one state
        regimes(:)=1;
        
    end
    
end
%-------------------------------------
if condforkst
    
    if ~isempty(switch_rule)
        
        error('conditional forecasting with switching rule not yet implemented')
        
    end
    
    myshocks=[];
    
    if options.burn
        
        error('conditional forecasting and burn-in not allowed')
        
    end
    
    model_is_linear=solve_order==1;
    
    if model_is_linear
        % apply standard tools
        standard_conditional_forecasting_tools();
        
    else
        % apply nonlinear tools
        if all_known_states
            % apply nonlinear tools
            nonlinear_conditional_forecasting_tools();
            
        else
            
            error('randomly jumping from one state to another not allowed for conditional forecasting')
            
        end
        
    end
    
else
    
    % condition on shocks only: The shock database may not have the
    % appropriate size if not all shocks are conditioned upon. This has to
    % be amended for.
    bigt=size(shocks,2);
    
    tmp=zeros(nx,bigt);
    
    tmp(y0.econd.pos,:)=shocks(:,:,1);
    
    % expand as necessary if the conditions do not span the whole
    % simulation period.
    tmp(:,bigt+1:span+options.k_future)=0;
    
    shocks__=condition_on_shocks_only(tmp);
    
    myshocks=shocks__(:,options.burn+1:end);
    
end

regimes=regimes(options.burn+1:end);

    function nonlinear_conditional_forecasting_tools()
        % here we do not substract the steady state as we do in the linear
        % case
        %------------------------------------------------------------------
        EndogenousConditions=recreate_conditions(y_conditions);
        
        ct=EndogenousConditions{1};
        
        horizon=options.k_future+1;
        
%         nap=horizon;
        shocks_=shocks;
        
        ncols_shocks=size(shocks_,2);
        
        ncols_ct=size(ct,2);
        
        if strcmpi(forecast_conditional_hypothesis,'nas')
            
            shocks_=shocks_(:,1:horizon);
            
            ct=ct(:,1:min(ncols_ct,horizon));
            
            ncols_shocks=size(shocks_,2);
            
            ncols_ct=size(ct,2);
            
        end
        
        needed_shocks=options.nsteps+options.k_future;
        
        missing_shocks=max(0,needed_shocks-ncols_shocks);%<--missing_shocks=max(0,ncsp-ncols_shocks);
        % set shocks beyond the shocks availability to zero.
        myshocks=[shocks_,zeros(nx,missing_shocks)];
        
        % the nan shocks are the ones to estimate
        %-----------------------------------------
        estim_shocks=isnan(myshocks);
        % Set the nan shocks to 0: a simple simulation with shocks will be
        % performed in the absence of further constraints. otherwise, the
        % shocks will be updated given 0 as initial conditions
        %----------------------------------------------------------
        myshocks(estim_shocks)=0;
        
        % trim the conditional information up to the number of steps
        %-----------------------------------------------------------
        endo_horizon=min(ncols_ct,options.nsteps);
        
        ct=ct(:,1:endo_horizon);
        
        % the non-nan elements are the targets to hit
        %--------------------------------------------
        targets=~isnan(ct);
        
        % restricted rows
        %----------------
        endo_restricted_rows=EndogenousConditions{end};
        % minimize the distance between the endogenous conditions and
        % the predictions by adjusting the unconstrained shocks
        %--------------------------------------------------------------
        x1=fsolve(@model_distance,myshocks(estim_shocks),...
            optimset('display','none','TolFun',sqrt(eps)));
%         x1=lsqnonlin(@model_distance,myshocks(estim_shocks),[],[],...
%             optimset('display','none','TolFun',sqrt(eps)));
        myshocks(estim_shocks)=x1;
        
        % rebuild the final sims
        %-----------------------
        condition_on_shocks_only(myshocks);
        
        function crit=model_distance(ex0)
            
            ex=myshocks;
            
            ex(estim_shocks)=ex0;
            
            condition_on_shocks_only(ex);
            
            yf=sims(endo_restricted_rows,:);
            
            crit=yf(targets)-ct(targets);
            
            bad=isinf(crit);
            
            crit(bad)=penalty*sign(crit(bad));
            
            bad=isnan(crit);
            
            crit(bad)=penalty;
            
        end
        
    end

    function standard_conditional_forecasting_tools()
        
        if options.simul_frwrd_back_shoot
            
            [cfkst,myshocks,PAI,retcode]=...
                utils.forecast.simul_forward_back_shooting(T,ss,y0,...
                state_vars_location,h,nx,Qfunc,regimes,options);
            
        else
            
            Tstar=T;
            
            if options.simul_anticipate_zero
                
                Tstar=NotAnt_T;
                
            end
            
            model=struct('T',{Tstar},'sstate',{ss},'state_cols',state_vars_location,...
                'Qfunc',Qfunc,'k',options.k_future,'nshocks',nx);
            
            opt=utils.miscellaneous.reselect_options(options,@utils.forecast.rscond.forecast);
            
            if h==1||~(isempty(regimes)||any(isnan(regimes)))
                
                [myshocks,regimes,PAI,retcode,cfkst,Qt]=utils.forecast.rscond.forecast(model,y0.y,...
                    y0.ycond,y0.econd,opt,regimes);
                
            else
                
                [myshocks,regimes,PAI,retcode,cfkst,Qt]=utils.forecast.rscond.loop_forecast(model,y0.y,...
                    y0.ycond,y0.econd,opt,regimes);
                
            end
            
        end
        
        % skip the initial conditions and add the mean
        if retcode
            
            sims=[];
            
        else
            % remove one period of history as it will be added back later
            % on
            sims=cfkst(:,2:end,:);
            
        end
        
    end

    function shocks=condition_on_shocks_only(shocks)
        
        if isfield(options,'occbin') && ~isempty(options.occbin)
            
            do_occbin_sim()
            
            return
            
        end
        
        Tstar=T;
        
        if options.simul_anticipate_zero
            
            Tstar=NotAnt_T;
            
        end
        
        y1=[];
        
        % shocks that are nan are shocks that are not conditioned on
        %------------------------------------------------------------
        shocks(isnan(shocks))=0;
        
        y00=y0;
        
        for t=1:span
            
            if ~retcode
                % process shocks
                %----------------
                shocks_t=shocks(:,t+(0:options.k_future));
                
                % compute transition matrix and switching probabilities
                %------------------------------------------------------
                [Q,retcode]=Qfunc(y00.y);
                
                if ~retcode
                    
                    if t==1
                        
                        Q=full(Q);
                        
                        Qt=Q(:,:,ones(1,options.nsteps));
                        
                    end
                    
                    if t>options.burn && t<span
                        
                        Qt(:,:,t-options.burn)=Q;
                        
                    end
                    
                    smallrt=regimes(t);
                    
                    if isnan(smallrt)
                        % the state is not known
                        if t==1
                            % draw from initial distribution
                        else
                            % draw conditional on yesterday's state
                            PAI=Q(regimes(t-1),:);
                            
                        end
                        
                    end
                    
                    % compute the forecast
                    %---------------------
                    state_list=1:h;
                    
                    iter=0;
                    
                    while (isempty(y1) && iter < h) && ~retcode
                        % draw state and compute forecast
                        %--------------------------------
                        one_step();
                        
                        iter=iter+1;
                        
                    end
                    
                    if isempty(y1)
                        
                        retcode=703;
                        
                    end
                    
                    if retcode
                        
                        return
                        
                    end
                    
                    if isnan(regimes(t))
                        
                        regimes(t)=smallrt;
                        
                    end
                    
                    if t>options.burn
                        
                        sims(:,t-options.burn)=y1.y;
                        
                    end
                    
                    y00=y1;
                    
                    y1=[];
                    
                end
                
            end
            
        end
        
        function one_step()
            
            ok = true;
            
            if ~isempty(switch_rule)
                
                if t==1 && ~iscell(switch_rule)
                    
                    switch_rule={switch_rule};
                    
                end
                % evaluate all regimes
                %---------------------
                ysr=y00(ones(1,h));
                
                for ireg=1:h
                    
                    ysr(ireg)=utils.forecast.one_step_fbs(Tstar(:,ireg),...
                        y00,ss{ireg},state_vars_location,...
                        options.simul_sig,shocks_t,options.simul_order);
                    
                end
                
                ysry=[ysr.y];
                
                ok=switch_rule{1}(ysry(inv_order_var,:),smallrt,...
                    regimes(1:t-1),...
                    sims(inv_order_var,1:t-1),switch_rule{2:end});
                
                % zero probability for offenders
                %-------------------------------
                PAI(~ok)=0;
                
            end
            
            cp=rebuild_cp();
            
            if retcode
                
                return
                
            end
            
            if isnan(smallrt)
                
                lucky=find(cp>rand,1,'first')-1;
                
                smallrt=state_list(lucky);
                
            end
            
            if ~isempty(switch_rule)
                
                y1=ysr(smallrt);
                
                return
                
            end
            
            myorder={Tstar(:,smallrt),y00,ss{smallrt},...
                state_vars_location,options.simul_sig,shocks_t,...
                options.simul_order};
            
            % use the solution prescribed by simul_anticipate_zero
            y1=utils.forecast.one_step_fbs(myorder{:});
                
            if ~isempty(options.complementarity) && ~options.complementarity(y1.y)
                
                if options.simul_honor_constraints_through_switch
                    
                    if ~isnan(regimes(t))
                        
                        warning(['forced to apply solution %0.0f ',...
                            'but cannot find a feasible path'],regimes(t)) 
                        
                        retcode=703; 
                                            
                    end
                    
                    ok = false;
                
                else
                    % use the normal solution
                    %-------------------------
                    [y1,~,retcode,shocks_t1]=utils.forecast.one_step_fbs(T(:,smallrt),y00,ss{smallrt},state_vars_location,...
                        options.simul_sig,shocks_t,options.simul_order,...
                        options.sep_compl,cond_shocks_id);
                    
                    if options.simul_anticipate_zero
                        % update contemporaneous shocks only since future
                        % anticipated shocks may change later because of
                        % unanticipated shocks... and we do not want to
                        % push zeros for the non conditional shocks.
                        shocks(:,t)=shocks_t1(:,1);
                    
                    else
                        
                        shocks(:,t+(0:options.k_future))=shocks_t1;
                    
                    end
                    
                end
                
            end
            
            if ~ok
                
                state_list(state_list==smallrt)=[];
                
                smallrt=nan;
                
                y1=[];
                
            end
            
            function cp=rebuild_cp()
                
                PAI00=PAI(state_list);
                
                PAI00=PAI00/sum(PAI00);
                
                if any(isnan(PAI00))
                    
                    retcode=703;
                    
                end
                
                cp=cumsum(PAI00);
                
                cp=[0,cp(:).'];
            
            end
            
        end
        
        function do_occbin_sim()
            
            [sims,regimes,retcode]=utils.forecast.simul_occbin(y0,T,ss,state_vars_location,...
                options,shocks);
            
            Q=full(Qfunc(y0.y));
            
            Qt=Q(:,:,ones(1,numel(regimes)-options.burn));
            
            sims=sims(:,options.burn+1:end);
        
        end
        
    end

    function NotAnt_T=no_anticipation_solution()
        
        NotAnt_T=T;
        
        if options.k_future==0
            
            return
            
        end
                
        nz=size(NotAnt_T{1,1},2);
        
        zproto=false(1,nz);
        
        offset=nsv+1+nx;
        
        bad=1:nx;
        
        for iplus=1:options.k_future
            
            bad_locs=offset+bad;
            
            zproto(bad_locs)=true;
            
            offset=offset+nx;
            
        end
        
        zkz=zproto;
        
        for io=1:solve_order
            
            for ireg_=1:h
                
                NotAnt_T{io,ireg_}(:,zkz)=0;
                
            end
            
            if io<solve_order
                
                zkz=logical(kron(zkz,zproto));
                
            end
            
        end
        
    end

end

function c=recreate_conditions(c,hard)

if nargin<2
    
    hard=true;
    
end

% remove the possibly singleton dimension
c=c(:,:);

% keep only the good rows
good=any(~isnan(c),2);

rest_id=find(good);

c=c(good,:);% <---c=permute(c(good,:),[2,1]);

if hard
    
    lb=c;
    
    ub=c;
    
else
    
    lb=-inf(size(c));
    
    ub=inf(size(c));
    
end

lbub=struct('LB',lb,'UB',ub);

OMG =[];

c={c,OMG,lbub,rest_id};

end