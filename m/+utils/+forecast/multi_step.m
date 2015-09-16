function [sims,states,retcode,Qt,myshocks]=multi_step(y0,ss,T,state_vars_location,options)
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

% Qt below may be used to get the time series of the transition matrix
endo_nbr=size(y0.y,1);
PAI=options.PAI;
Qfunc=options.Qfunc;
states=options.states(:);
shocks=options.shocks;
cond_shocks_id=[];
if options.k_future
    if isfield(options,'shock_structure')
        cond_shocks_id=any(options.shock_structure,2);
    else
        cond_shocks_id=any(isnan(shocks(:,2:end)),2);
    end
end

simul_update_shocks_handle=options.simul_update_shocks_handle;
simul_do_update_shocks=options.simul_do_update_shocks;
forecast_conditional_hypothesis=options.forecast_conditional_hypothesis;
options=rmfield(options,{'states','shocks','PAI','Qfunc','y'});
y_conditions=[];
if size(y0.y,3)>1
    y_conditions=y0.y(:,:,2:end);
    y0.y=y0.y(:,:,1);
end
h=size(T,2);
Qt=[];
retcode=0;
penalty=1e+6;

do_Qt=nargout>3;
span=options.nsteps+options.burn;
sims=nan(endo_nbr,options.nsteps);
nsv=numel(state_vars_location);
nx=(size(T{1,1},2)-nsv-1)/(options.k_future+1); %<---size(shocks,1);

condforkst=~isempty(y_conditions);
%-------------------------------------
all_known_states=~any(isnan(states));
unique_state=numel(PAI)==1||(all_known_states && all(states==states(1)));
% model_is_linear=size(T,1)==1;
if unique_state
    if isnan(states(1))
        % no states were provided. Uniqueness in this case implies that
        % we have only one state
        states(:)=1;
    end
end
%-------------------------------------
if condforkst
    myshocks=[];
    if options.burn
        error('conditional forecasting and burn-in not allowed')
    end
    model_is_linear=size(T,1)==1;
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
    % condition on shocks only
    shocks__=condition_on_shocks_only(shocks);
    myshocks=shocks__(:,options.burn+1:end);
 end

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
        model=struct('T',{T},'sstate',{ss},'state_cols',state_vars_location,...
            'Qfunc',Qfunc,'k',options.k_future,'nshocks',nx);
        opt=utils.miscellaneous.reselect_options(options,@utils.forecast.rscond.forecast);
        if h==1||~(isempty(states)||any(isnan(states)))
            [myshocks,states,PAI,retcode,cfkst]=utils.forecast.rscond.forecast(model,y0.y,...
                y0.ycond,y0.econd,opt,states);
        else
            [myshocks,states,PAI,retcode,cfkst]=utils.forecast.rscond.loop_forecast(model,y0.y,...
                y0.ycond,y0.econd,opt,states);
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
                if ~isempty(simul_update_shocks_handle) && simul_do_update_shocks
                    shocks_t=simul_update_shocks_handle(shocks_t,y00.y);
                end
                if options.simul_anticipate_zero
                    shocks_t(:,2:end)=0;
                end
                % compute transition matrix and switching probabilities
                %------------------------------------------------------
                [Q,retcode]=Qfunc(y00.y);
                if ~retcode
                    rt=states(t);
                    if isnan(rt)
                        % the state is not known
                        if t==1
                            % draw from initial distribution
                        else
                            % draw conditional on yesterday's state
                            PAI=Q(states(t-1),:);
                        end
                        if do_Qt
                            if t==1
                                Qt=Q(:,:,ones(1,options.nsteps));
                            end
                            if t>options.burn && t<span
                                Qt(:,:,t-options.burn)=Q;
                            end
                        end
                    end
                    
                    % compute the forecast
                    %---------------------
                    state_list=1:h;
                    iter=0;
                    while isempty(y1) && iter < h
                        % draw state and compute forecast
                        %--------------------------------
                        one_step();
                        iter=iter+1;
                    end
                    if isempty(y1)
                        error('I could not find a feasible path')
                    end
                    if isnan(states(t))
                        states(t)=rt;
                    end
                    
                    if t>options.burn
                        sims(:,t-options.burn)=y1.y;
                    end
                    y00=y1;
                    y1=[];
                end
            end
        end
        states=states(options.burn+1:end);
        function one_step()
            cp=rebuild_cp();
            if isnan(rt)
                lucky=find(cp>rand,1,'first')-1;
                rt=state_list(lucky);
            end
            y1=utils.forecast.one_step_fbs(T(:,rt),y00,ss{rt},state_vars_location,...
                options.simul_sig,shocks_t,options.simul_order);
            if ~isempty(options.complementarity) && ~options.complementarity(y1.y)
                if options.simul_honor_constraints_through_switch
                    if ~isnan(states(t))
                        error(sprintf('forced to apply solution %0.0f but cannot find a feasible path',states(t))) %#ok<SPERR>
                    end
                    state_list(state_list==rt)=[];
                    rt=nan;
                    y1=[];
                else
                    [y1,~,retcode,shocks_t1]=utils.forecast.one_step_fbs(T(:,rt),y00,ss{rt},state_vars_location,...
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
            function cp=rebuild_cp()
                PAI00=PAI(state_list);
                PAI00=PAI00/sum(PAI00);
                if any(isnan(PAI00))
                    error('I could not find a feasible path')
                end
                cp=cumsum(PAI00);
                cp=[0,cp(:).'];
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