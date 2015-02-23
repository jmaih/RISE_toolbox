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
    if unique_state
        if model_is_linear
            % apply standard tools
            standard_conditional_forecasting_tools();
        else
            % apply nonlinear tools
            nonlinear_conditional_forecasting_tools();
        end
    else
        if all_known_states
            % apply nonlinear tools
            nonlinear_conditional_forecasting_tools();
        else
            error('randomly jumping from one state to another not allowed for conditional forecasting')
        end
    end
else
    % condition on shocks only
    condition_on_shocks_only(shocks);
    myshocks=shocks;
end

    function nonlinear_conditional_forecasting_tools()
        % here we do not substract the steady state as we do in the linear
        % case
        %------------------------------------------------------------------
        EndogenousConditions=recreate_conditions(y_conditions);
        ct=transpose(EndogenousConditions{1});
%         if any(any(shocks))
%             ShocksConditions=recreate_conditions(shocks);
%         else
%             ShocksConditions=[];
%         end
        horizon=options.k_future+1;
        nap=horizon;
        ncp=size(ct,2);
        ncsp = utils.forecast.conditional.number_of_conditioning_shocks_periods(...
            forecast_conditional_hypothesis,ncp,nap);
        missing_shocks=max(0,ncsp-size(shocks,2));
        myshocks=[shocks,nan(nx,missing_shocks)];
        myshocks=myshocks(:,1:ncsp);
        % set all zero shocks to nans : we probably don't want to condition
        % on zero shocks, this is most likely the result of initialization
        %------------------------------------------------------------------
        myshocks(myshocks==0)=nan;
        % the nan shocks are the ones to estimate
        estim_shocks=isnan(myshocks);
        % Set the nan shocks to 0: a simple simulation with shocks will be
        % performed in the absence of further constraints. otherwise, the
        % shocks will be updated given 0 as initial conditions
        %----------------------------------------------------------
        myshocks(estim_shocks)=0;
        
        % trim the conditional information up to the number of steps
        %-----------------------------------------------------------
        endo_horizon=min(ncp,options.nsteps);
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
            optimset('display','iter','TolFun',sqrt(eps)));
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
        rt=states(1);
        % re-build the square matrix
        %----------------------------
        n=size(T{rt},1);
        H=zeros(n);
        H(:,state_vars_location)=T{rt}(:,1:nsv);
        risk=T{rt}(:,nsv+1);
        G=reshape(T{rt}(:,nsv+2:end),[n,nx,options.k_future+1]);
        ss_=ss{rt};
        if any(risk)
            ss_=ss_+(eye(n)-H)\risk;
        end
        Y0=bsxfun(@minus,y0.y,ss_);
        NumberOfSimulations=0;
        EndogenousConditions=recreate_conditions(bsxfun(@minus,y_conditions,ss_));
        if ~isempty(shocks)
            ShocksConditions=recreate_conditions(shocks);
        else
            ShocksConditions=[];
        end
        [~,CYfMean,myshocks]=utils.forecast.conditional.forecast_engine(...
            Y0,H,G,EndogenousConditions,ShocksConditions,options.nsteps,...
            NumberOfSimulations,forecast_conditional_hypothesis);
        % skip the initial conditions and add the mean
        sims=bsxfun(@plus,CYfMean(:,2:end),ss_);
    end

    function condition_on_shocks_only(shocks)
        y1=[];
        for t=1:span
            if ~retcode
                % process shocks
                %----------------
                shocks_t=shocks(:,t+(0:options.k_future));
                if ~isempty(simul_update_shocks_handle) && simul_do_update_shocks
                    shocks_t=simul_update_shocks_handle(shocks_t,y0.y);
                end
                
                % compute transition matrix and switching probabilities
                %------------------------------------------------------
                [Q,retcode]=Qfunc(y0.y);
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
                    y0=y1;
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
            y1=utils.forecast.one_step(T(:,rt),y0,ss{rt},state_vars_location,...
                options.simul_sig,shocks_t,options.simul_order);
            if ~options.complementarity(y1.y)
                if ~isnan(states(t))
                    error(sprintf('forced to apply solution %0.0f but cannot find a feasible path',states(t))) %#ok<SPERR>
                end
                state_list(state_list==rt)=[];
                rt=nan;
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

function c=recreate_conditions(c)
% remove the possibly singleton dimension
c=c(:,:);
% keep only the good rows
good=any(~isnan(c),2);
rest_id=find(good);
c=c(good,:);% <---c=permute(c(good,:),[2,1]);
lb=-inf(size(c));
ub=inf(size(c));
lbub=struct('LB',lb,'UB',ub);
OMG =[];
c={c,OMG,lbub,rest_id};
end