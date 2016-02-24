function [sims,regimes,retcode,Qt,myshocks]=multi_step(y0,ss,T,state_vars_location,options)
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

% use first page only as it is the most relevant
regimes=vec(y0.rcond.data(:,:,1));

shocks=y0.econd.data;
cond_shocks_id=[];
if options.k_future
    if isfield(options,'shock_structure')
        cond_shocks_id=any(options.shock_structure,2);
    else
        cond_shocks_id=any(isnan(shocks(:,2:end,1)),2);
    end
end

simul_update_shocks_handle=options.simul_update_shocks_handle;
simul_do_update_shocks=options.simul_do_update_shocks;
forecast_conditional_hypothesis=options.forecast_conditional_hypothesis;
options=rmfield(options,{'PAI','Qfunc','y'});

[solve_order,h]=size(T);
Qt=[];
retcode=0;
penalty=1e+6;

do_Qt=nargout>3;
span=options.nsteps+options.burn;
sims=nan(endo_nbr,options.nsteps);
nsv=numel(state_vars_location);
nx=(size(T{1,1},2)-nsv-1)/(options.k_future+1); %<---size(shocks,1);

condforkst=~isempty(y0.ycond.data)||...
    (...
    ~isempty(shocks) && (any(isnan(shocks(:)))||any(isinf(shocks(:)))||...
    max(max(abs(shocks(:,:,1)-shocks(:,:,2))>1e-19))||...
    max(max(abs(shocks(:,:,1)-shocks(:,:,3))>1e-19))||...
    max(max(abs(shocks(:,:,2)-shocks(:,:,3))>1e-19)))...
    )||...
    options.forecast_conditional_sampling_ndraws>1;

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
        Tstar=T;
        if options.simul_anticipate_zero
            Tstar=NotAnt_T;
        end
        model=struct('T',{Tstar},'sstate',{ss},'state_cols',state_vars_location,...
            'Qfunc',Qfunc,'k',options.k_future,'nshocks',nx);
        opt=utils.miscellaneous.reselect_options(options,@utils.forecast.rscond.forecast);
        if h==1||~(isempty(regimes)||any(isnan(regimes)))
            [myshocks,regimes,PAI,retcode,cfkst]=utils.forecast.rscond.forecast(model,y0.y,...
                y0.ycond,y0.econd,opt,regimes);
        else
            [myshocks,regimes,PAI,retcode,cfkst]=utils.forecast.rscond.loop_forecast(model,y0.y,...
                y0.ycond,y0.econd,opt,regimes);
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
                if ~isempty(simul_update_shocks_handle) && simul_do_update_shocks
                    shocks_t=simul_update_shocks_handle(shocks_t,y00.y);
                end
                % compute transition matrix and switching probabilities
                %------------------------------------------------------
                [Q,retcode]=Qfunc(y00.y);
                if ~retcode
                    rt=regimes(t);
                    if isnan(rt)
                        % the state is not known
                        if t==1
                            % draw from initial distribution
                        else
                            % draw conditional on yesterday's state
                            PAI=Q(regimes(t-1),:);
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
                    if isnan(regimes(t))
                        regimes(t)=rt;
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
            cp=rebuild_cp();
            if isnan(rt)
                lucky=find(cp>rand,1,'first')-1;
                rt=state_list(lucky);
            end
            % use the solution prescribed by simul_anticipate_zero
            y1=utils.forecast.one_step_fbs(Tstar(:,rt),y00,ss{rt},state_vars_location,...
                options.simul_sig,shocks_t,options.simul_order);
            if ~isempty(options.complementarity) && ~options.complementarity(y1.y)
                if options.simul_honor_constraints_through_switch
                    if ~isnan(regimes(t))
                        error(sprintf('forced to apply solution %0.0f but cannot find a feasible path',regimes(t))) %#ok<SPERR>
                    end
                    state_list(state_list==rt)=[];
                    rt=nan;
                    y1=[];
                else
                    % use the normal solution
                    %-------------------------
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
        function do_occbin_sim()
            [sims,regimes,retcode]=simul_occbin(y0,T,ss,state_vars_location,...
                options,shocks);
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



function [sim1,regimes,retcode]=simul_occbin(y0,T,ss,state_vars_location,...
    options,shocks)

use_pinv=false;
is_accelerated=true;

% find the first-order approximation of the system
%--------------------------------------------------

A0=options.occbin.A0;
Aplus=A0;
h=size(A0,3);
for ireg=1:h
    Aplus(:,:,ireg)=options.occbin.Gplus01(:,:,ireg,ireg);
end
Aminus=options.occbin.Aminus;
B=options.occbin.B;
c=options.occbin.c;

ref_state=options.solve_occbin;
other_state=setdiff(1:h,ref_state);
if numel(other_state)>1
    error('occbin can only handle one constraint')
end
nstate_vars=numel(state_vars_location);

endo_nbr=size(y0.y,1);
exo_nbr=size(B,2);

span=options.nsteps+options.burn;

compl=@(x)isempty(options.complementarity)||options.complementarity(x);

% find the solution for the reference regime
%--------------------------------------------
[H,G,k]=reference_solution();

retcode=0;
y0.econd.data=shocks(:,:,ones(3,1));

last_step=0;
bigt=0;
Ht0=[];
Gt0=[];
kt0=[];

% initialize with nans in order to see if/when simulation breaks down
sim1=nan(endo_nbr,span);
regimes=ref_state*ones(span,1);
do_one_occbin_path()

    function [H,G,k]=reference_solution()
        H=zeros(endo_nbr);
        H(:,state_vars_location)=T{ref_state}(:,1:nstate_vars);
        k=full(T{ref_state}(:,nstate_vars+1));
        G=full(T{ref_state}(:,nstate_vars+2:end));
    end

    function do_one_occbin_path()
        % try and forecast from the ref
        %------------------------------
        Tfunc=@(t)T{ref_state};
        is_viol=forecaster(y0.y,Tfunc);
        if ~is_viol
            return
        end
        bigt=1;
        last_step=1;
        % accelerate things by preparing the longest possible stretch
        %-------------------------------------------------------------
        if is_accelerated
            [Ht0,Gt0,kt0]=resolve_occbin(span-1);
        end
        
        % make assumption about the duration of the violation and retry
        %---------------------------------------------------------------
        if last_step<=1
            y00=y0.y;
        else
            y00=sim1(:,last_step-1);
        end
        while is_viol && last_step-1+bigt<span
            [Ht,Gt,kt]=load_solutions(bigt);
            if retcode
                % we exit as quickly as possible, it is not possible to
                % tame the system.
                bigt=inf;
            else
                Tfunc=@(t)[Ht(:,state_vars_location,t),kt(:,t),Gt(:,:,t)];
                is_viol=forecaster(y00,Tfunc,last_step);
                if is_viol
                    bigt=bigt+1;
                    if last_step<=1
                        y00=y0.y;
                    else
                        y00=sim1(:,last_step-1);
                    end
                end
            end
        end
        if is_viol
            retcode=701;
        end
    end

    function [Ht,Gt,kt]=load_solutions(bigt)
        if is_accelerated
            last_step_=max(1,last_step);
            wing=last_step_:span;
            nwing=numel(wing);
            stretch=span-bigt:span-1;
            Ht=H(:,:,ones(nwing,1));Ht(:,:,1:bigt)=Ht0(:,:,stretch);
            Gt=G(:,:,ones(nwing,1));Gt(:,:,1:bigt)=Gt0(:,:,stretch);
            kt=k(:,ones(nwing,1));kt(:,1:bigt)=kt0(:,stretch);
            regimes(last_step_:last_step_+bigt-1)=other_state;
            regimes(last_step_+bigt:end)=ref_state;
        else
            [Ht,Gt,kt]=resolve_occbin(bigt);
        end
    end

    function is_viol=forecaster(y00,Tfunc,istart)
        if nargin<3
            istart=1;
        end
        iter=0;
        for istep=istart:span
            iter=iter+1;
            y00_s=y00-ss{1}; % steady state is the same???
            state=[y00_s(state_vars_location);1;shocks(:,istep)];
            sim1(:,istep)=ss{1}+Tfunc(istep-istart+1)*state;
            y00=sim1(:,istep);
            if any(~isfinite(y00))||~isreal(y00)
                retcode=701;
            end
            is_viol=~compl(y00);
            if is_viol||retcode
                break
            end
        end
        % if the period after the guess does not violate, then the guess
        % was correct. Update the last step
        %----------------------------------------------------------------
        if ~is_viol||(iter>bigt && compl(sim1(:,istart+bigt)))
            last_step=istep;
            bigt=0;
        end
    end

    function [Ht,Gt,kt]=resolve_occbin(bigt)
        new_span=span-last_step+1;
        Ht=zeros(endo_nbr,endo_nbr,new_span);
        Gt=zeros(endo_nbr,exo_nbr,new_span);
        kt=zeros(endo_nbr,new_span);
        
        t=bigt+1:new_span;
        Ht(:,:,t)=H(:,:,ones(1,new_span-bigt));
        Gt(:,:,t)=G(:,:,ones(1,new_span-bigt));
        kt(:,t)=k(:,ones(1,new_span-bigt));
        regimes(t+last_step-1)=ref_state;
        
        maxed_out=bigt>=new_span;
        if maxed_out
            retcode=701;
        else
            for t=bigt:-1:1
                if use_pinv
                    AHAi=-pinv(Aplus(:,:,other_state)*Ht(:,:,t+1)+...
                        A0(:,:,other_state));
                else
                    AHAi=-(Aplus(:,:,other_state)*Ht(:,:,t+1)+...
                        A0(:,:,other_state))\eye(endo_nbr);
                end
                Ht(:,:,t)=AHAi*Aminus(:,:,other_state);
                Gt(:,:,t)=AHAi*B(:,:,other_state);
                kt(:,t)=AHAi*(c(:,other_state)+Aplus(:,:,other_state)*kt(:,t+1));
                regimes(t+last_step-1)=other_state;
            end
        end
    end
end