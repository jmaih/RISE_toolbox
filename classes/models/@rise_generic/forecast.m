function cond_fkst_db=forecast(obj,varargin)% historical_db
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% hist_db,conditions,ShockProperties,uncertainty,Nsteps,Nsim,Hypothesis,method)
% HistDB,CondDB,OMGy,OMGx
% 'forecast_conditional_info' to be structure with fields:
% lower_bound,upper_bound,central_tendency, all of them being
% ts

if isempty(obj)
    if nargout>1
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    end
    cond_fkst_db=struct('forecast_conditional_hypothesis',0,...
        'forecast_conditional_info',[],...
        'forecast_start_date','',...
        'forecast_nsteps',12);
    %         'truncated_normal_simulation_method','ghk'
    % alternative is 'gibbs'% should be moved elsewhere
    
    % OMGy,OMGx
    % alternatives for forecast_conditional_hypothesis are:
    % 'ncp'(or 'iris'): Number of conditioning periods
    % 'nas'(or 'iwb'): Number of anticipated steps
    % 0 (or 'jma'): No hypothesis
    
    return
end
%% pass the options
obj=set(obj,varargin{:});

%% pass some options to the simulate method
if isempty(obj.options.forecast_start_date)
    if isempty(obj.options.estim_end_date)
        error('cannot find a date to start the forecast')
    end
    obj=set(obj,'forecast_start_date',serial2date(date2serial(obj.options.estim_end_date)+1));
end

end_date='';
if ~isempty(obj.options.forecast_start_date)
    end_date=serial2date(date2serial(obj.options.forecast_start_date)-1);
end
fkstdata='';
if ~isempty(obj.options.data)
    fkstdata=obj.options.data;
end
obj=set(obj,...
    'simul_periods',obj.options.forecast_nsteps,...
    'simul_history_end_date',end_date,...
    'simul_historical_data',fkstdata);

cond_fkst_db=simulate(obj);

return
%% re-collect the options
if obj.markov_chains.regimes_number==1
    obj.options.forecast_regime=1;
end

% truncated_normal_simulation_method=obj.options.truncated_normal_simulation_method;
forecast_conditional_hypothesis=obj.options.forecast_conditional_hypothesis;

%% solution and system
%---------------------
[obj,retcode]=solve(obj);
if retcode
    error([mfilename,':: no solution: the model cannot be simulated'])
end

%% initial conditions
% note that initial conditions for the probabilities should also be given
% by the smoother in case of forecasting. This is something that should
% implemented later on
%--------------------------------------------------------------------------
Initcond=set_simulation_initial_conditions(obj);
y0=Initcond.y.y;
PAI=Initcond.PAI;
% exogenous_observed=Initcond.exogenous_observed;
simul_burn=Initcond.burn;
endo_nbr=obj.endogenous.number(end);
Qfunc=Initcond.Qfunc;
y0cols=size(y0,2);

% if isempty(obj.options.simul_histval)
%     error('historical values not provided')
% end
%% conditional information

conditional_info=obj.options.forecast_conditional_info;
forecast_start_date=serial2date(1+date2serial(Initcond.simul_history_end_date));
conditions=conditional_forecasting.parse_information(conditional_info,...
    obj.endogenous.name,obj.exogenous.name,date2serial(forecast_start_date));

%%
h=obj.markov_chains.regimes_number;
% number of forecasting steps
nsteps=obj.options.simul_periods;
simul_regime=obj.options.simul_regime;
if ~isempty(simul_regime) && ~all(ismember(simul_regime,1:h))
    error([mfilename,':: all elements in the history of the markov states ',...
        'must be between 1 and ',int2str(h)])
end

% generate the state vector, which will be used also under optimization
%----------------------------------------------------------------------
State=generic_tools.set_simulation_regimes(simul_regime,nsteps,simul_burn);

nexo=sum(obj.exogenous.number);
if isempty(conditions)
    % generate zero shocks: it is enough to generate one observation
    %---------------------------------------------------------------
    if obj.options.simul_no_shocks
        shocks=[];
    else
        shocks=zeros(nexo,1);
    end
else
    %     if ((isa(obj,'dsge') && obj.options.solve_order ==1)||isa(obj,'svar')) && ...
    %             h==1||(~isempty(simul_regime) && all(simul_regime)==simul_regime(1))
    %         % apply the old algorithm and return immediately
    %         %-----------------------------------------------
    %         keyboard
    %     else
    ct=conditions.central_tendency;
    % separate the constrained shocks from the unconstrained ones
    %------------------------------------------------------------
    n_conditioning_periods=size(ct,1);
    shocks=nan(n_conditioning_periods,nexo);
    locs=conditions.restrictions_id(~conditions.endogenous_vars);
    shocks(:,locs)=ct(:,~conditions.endogenous_vars);
    % trim the conditional information: future shocks can exceed nsteps
    % only in the dsge case and we will have nsteps+horizon-1 as the
    % maximum horizon for the shocks
    %------------------------------------------------------------------
    horizon=1;
    if isa(obj,'dsge')
        horizon=max(obj.exogenous.shock_horizon);
    end
    shocks_length=nsteps+horizon-1;
    if isa(obj,'dsge') && horizon>1
        switch lower(forecast_conditional_hypothesis)
            case {'ncp','iwb'}%NumberOfCondShocksPeriods=ncp;
                shocks_length=min(shocks_length,n_conditioning_periods);
            case {'nas','iris'}
                shocks_length=horizon;
                if ~isequal(shocks_length,n_conditioning_periods)
                    error([mfilename,':: for the NAS or IRIS assumption, you need # anticipated steps = # conditioning periods'])
                end
            case {0,'jma'}
                % done
            otherwise
                error([mfilename,':: Unknown option for the anticipation hypothesis'])
        end
    end
    missing=max(0,shocks_length-n_conditioning_periods);
    shocks=transpose([shocks;nan(missing,nexo)]);
    shocks=shocks(:,1:shocks_length);
    % the nan shocks are the ones to estimate
    estim_shocks=isnan(shocks);
    % Set the nan shocks to 0: a simple simulation with shocks will be
    % performed in the absence of further constraints. otherwise, the
    % shocks will be updated given 0 as initial conditions
    %----------------------------------------------------------
    shocks(estim_shocks)=0;
    if ~isempty(ct)
        % trim the conditional information up to the number of steps
        %-----------------------------------------------------------
        ct(:,~conditions.endogenous_vars)=[];
        endo_horizon=min(n_conditioning_periods,nsteps);
        ct=ct(1:endo_horizon,:);
        ct=transpose(ct);
        % the non-nan elements are the targets to hit
        %--------------------------------------------
        targets=~isnan(ct);
        % restricted rows
        %----------------
        endo_restricted_rows=conditions.restrictions_id(conditions.endogenous_vars);
        % minimize the distance between the endogenous conditions and
        % the predictions by adjusting the unconstrained shocks
        %--------------------------------------------------------------
        % x1=fminunc(@norm_model_distance,shocks(estim_shocks),optimset('display','iter','TolFun',sqrt(eps)));
        x1=fsolve(@model_distance,shocks(estim_shocks),...
            optimset('display','none','TolFun',sqrt(eps)));
        shocks(estim_shocks)=x1;
    end
    %     end
end

%% generate forecasts conditional on shocks
%-----------------------------------------
y=compute_forecasts_conditional_on_shocks(shocks);

%% store the simulations in a database
%------------------------------------
hist_start_date=serial2date(date2serial(Initcond.simul_history_end_date)-y0cols+1);
cond_fkst_db=ts(hist_start_date,y',obj.endogenous.name);
cond_fkst_db=pages2struct(cond_fkst_db);
%% also store the shocks if available
%------------------------------------
if ~isempty(shocks)
    shock_db=ts(forecast_start_date,shocks',obj.exogenous.name);
    cond_fkst_db=utils.miscellaneous.mergestructures(cond_fkst_db,pages2struct(shock_db));
end
% obj.forecasts=ts(obj.options.data.start,permute(y,[]),char(obj.endogenous.name));

% I need to select the conditional information. I still need to make up my
% mind on how to deal with the dichotomy of real_time vs soft conditions

% another issue

% and then splice the historical database and the forecast database
%

    function y=compute_forecasts_conditional_on_shocks(shocks)
        y0_linear=[];
        y0_=y0;
        for t=1:nsteps
            % draw a state
            %-------------
            [State(t+1),Q0,PAI,retcode]=generic_tools.choose_state(State(t+1),Qfunc,PAI,y0_(1:endo_nbr,end));
            if retcode
                return
            end
            % Simulate one step: since with have shocks, zt is irrelevant
            %------------------------------------------------------------
            zt=[];
            [y1,y0_linear]=simulation_engine(obj,y0_,y0_linear,zt,State(t+1),shocks(:,t:end));
            if t==1
                y=[y0_(1:endo_nbr,:),zeros(endo_nbr,nsteps)];
            end
            y(:,t+y0cols)=y1(1:endo_nbr,end);
            y0_=y1;
        end
    end

%     function crit=norm_model_distance(ex0)
%         crit=model_distance(ex0);
%         crit=100*norm(crit);
%     end

    function crit=model_distance(ex0)
        ex=shocks;
        ex(estim_shocks)=ex0;
        y=compute_forecasts_conditional_on_shocks(ex);
        yf=y(endo_restricted_rows,y0cols+1:end);
        crit=yf(targets)-ct(targets);
    end

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

% function [rstar,rlow,rhigh,OMG]=ConditionalDistributionOfShocks(MU,OMG,LB,UB,DTy)
% rstar=MU-DTy;
% rlow=LB-DTy;
% rhigh=UB-DTy;
%
% end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

% function xx=TruncatedMultivariateNormalMean(m,lb,ub,CholCov)
% npar=size(CholCov,1);
% xx=nan(npar,1);
% BOUNDS=[lb-m,ub-m];
% for i=1:npar
%     tmp=(BOUNDS(i,:)-CholCov(i,1:i-1)*xx(1:i-1,1))/CholCov(i,i);
%     xx(i)=TruncatedNormalMean(tmp(1),tmp(2));
% end
% xx=m+CholCov*xx;
% end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

% function y=TruncatedNormalMean(lb,ub)
% PHIl = NormalCumulativeDistribution(lb);
% PHIr = NormalCumulativeDistribution(ub);
% crit=1e-8;
%
% finite_lb_and_ub=isfinite(lb) & isfinite(ub);
% same=finite_lb_and_ub &  abs(ub-lb)<crit;
% tails=abs(PHIr-PHIl)<crit & ~ same;
% good_tails=tails & finite_lb_and_ub;
% bad_tails=tails & ~finite_lb_and_ub;
% others=~same & ~tails; clear tails
%
% if same % same, no problems
%     y=ub;
% elseif good_tails% assume a uniform distribution in the tails
%     y=.5*(ub+lb);
% elseif others% normal distribution for nice ones
%     y=-(NormalDensity(ub)-NormalDensity(lb))/(PHIr-PHIl);
% elseif bad_tails% Nasty ones
%     y=0;
% end
% % if abs(lb-ub)<1e-9
% %     y=ub;
% % else
% %     y=-(NormalDensity(ub)-NormalDensity(lb))/...
% %         (NormalCumulativeDistribution(ub)-NormalCumulativeDistribution(lb));
% % end
% end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

% function cdf=NormalCumulativeDistribution(x)
% cdf=.5*(1+erf(x/sqrt(2)));
% end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

% function pdf=NormalDensity(x)
% pdf=1/sqrt(2*pi)*exp(-.5*x^2);
% end


% forecast_regime,cond_data,cond_var_id,SS
% extract the steady state from the SS_and_BGP
% SS=obj.steady_state_and_balanced_growth_path(1:obj.endogenous.number(2),:);
% [cond_fkst,cond_fkst_mean,uncond_fkst,junk,junk]=...
%     forecast_engine(y0(:,regime_to_forecast),obj.T(:,:,regime_to_forecast),...
%     obj.R(:,:,:,regime_to_forecast),SS(:,regime_to_forecast));% compute the conditional forecasts

%     function [CYf,CYfMean,Yf,Emean,E]=forecast_engine(Y0,H,G,SS)
%         % Description:
%         % Usage:
%         % Arguments:
%         % Value:
%         % Author(s): Junior Maih (junior.maih@norges-bank.no)
%         % See Also:
%         % Examples:
%         NumberOfSimulations=1; % THIS IS TO REMAIN EXACTLY LIKE THIS !!!!
%
%         %% demean the initial condition.
%         yy0=Y0-SS;
%         %% demean the endogenous conditions
%         if ~isempty(EndogenousConditions)
%             rest_id	=EndogenousConditions{4};
%             SS_restr=transpose(SS(rest_id));
%             EndogenousConditions{1}=bsxfun(@minus,EndogenousConditions{1},SS_restr);% [ncp,ncv]=size(MU);
%             EndogenousConditions{3}.LB=bsxfun(@minus,EndogenousConditions{3}.LB,SS_restr);
%             EndogenousConditions{3}.UB=bsxfun(@minus,EndogenousConditions{3}.UB,SS_restr);
%         end
%         %%
%         [PbParams,M1,M2,RM2i,OMG,DYbar,junk,junk,Yf,MU,LB,UB]=...
%             ConditionalProjectionSubEngine(forecast_conditional_hypothesis,nsteps,H,G,yy0,...
%             EndogenousConditions,ShocksConditions);
%
%         [rstar,rlow,rhigh,OMG]=ConditionalDistributionOfShocks(MU,OMG,LB,UB,DYbar);
%
%         % check whether there are more holes
%         if any(isnan(rlow))||any(isnan(rhigh))
%             error([mfilename,':: LB and UB cannot have holes at other places than the holes in the Central tendency'])
%         end
%         if any(rlow>rhigh)
%             error([mfilename,':: some elements in LB greater than their counterparts in UB'])
%         end
%
%         CSOMG=transpose(chol(OMG));
%         % Compute the theoretical mean (doesn't change with the tightening in a symmetric
%         % distribution)
%         rmean=TruncatedMultivariateNormalMean(rstar,rlow,rhigh,CSOMG);
%         gam2mean=RM2i*rmean;
%         epsilon=M2*gam2mean;
%
%         ShockHorizon=PbParams.NumberOfAnticipatedPeriods-1+PbParams.maxiter;
%         Emean=InitializeArray(@zeros,PbParams.exo_nbr,ShockHorizon);
%         Emean(:,1:PbParams.NumberOfCondShocksPeriods)=reshape(epsilon,PbParams.exo_nbr,[]);
%         % Compute conditional forecasts
%         CYfMean=InitializeArray(@zeros,PbParams.endo_nbr,PbParams.maxiter,1,yy0);
%         CYfMean=ComputeForecasts(CYfMean,H,G,Emean,nsteps,PbParams.NumberOfAnticipatedPeriods);
%
%         seed=[];
%         r=TruncatedMultivariateNormalRnd(rstar,OMG,rlow,rhigh,...
%             NumberOfSimulations,truncated_normal_simulation_method,seed);
%         % draw gam2
%         gam2=RM2i*r;
%         clear r
%         % draw gam1
%         gam1=randn(PbParams.kk-PbParams.qq,NumberOfSimulations);
%         % form epsilon
%         epsilon=M1*gam1+M2*gam2;
%         clear gam1 gam2
%
%         % Build E
%         E=InitializeArray(@randn,PbParams.exo_nbr,ShockHorizon,NumberOfSimulations);
%         if ~isequal(forecast_conditional_hypothesis,0)
%             E=0*E;
%         end
%         E(:,1:PbParams.NumberOfCondShocksPeriods,:)=reshape(epsilon,[PbParams.exo_nbr,PbParams.NumberOfCondShocksPeriods,NumberOfSimulations]);
%         clear epsilon
%
%         CYf=InitializeArray(@zeros,PbParams.endo_nbr,PbParams.maxiter,NumberOfSimulations,yy0);
%         for isim=1:NumberOfSimulations
%             % Compute conditional forecasts
%             CYf(:,:,isim)=ComputeForecasts(CYf(:,:,isim),H,G,E(:,:,isim),nsteps,PbParams.NumberOfAnticipatedPeriods);
%         end
%
%         % Reset the forecasts to their normal length and re-mean them
%         CYf=bsxfun(@plus,CYf(:,1:nsteps+1,:),SS);
%         CYfMean=bsxfun(@plus,CYfMean(:,1:nsteps+1),SS);
%         Yf=bsxfun(@plus,Yf(:,1:nsteps+1),SS);
%         % and also the shocks (retain only the ones that matter)
%         ShockHorizonFinal=PbParams.NumberOfAnticipatedPeriods-1+nsteps;
%         E=E(:,1:ShockHorizonFinal,:);
%         Emean=Emean(:,1:ShockHorizonFinal);
%     end