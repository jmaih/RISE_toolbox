function [cond_fkst_db,cond_fkst_mean_db,uncond_fkst_db]=forecast(obj,varargin)% historical_db
% hist_db,conditions,ShockProperties,uncertainty,Nsteps,Nsim,Hypothesis,method)
% HistDB,CondDB,OMGy,OMGx
% 'forecast_conditional_info' to be structure with fields:
% lower_bound,upper_bound,central_tendency, all of them being
% rise_time_series  

if isempty(obj)
    if nargout>1
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    end
    cond_fkst_db=struct('forecast_conditional_hypothesis',0,...
        'forecast_paths',1000,...
        'forecast_param_uncert',true,...
        'forecast_shock_uncert',true,...
        'forecast_state_uncert',true,...
        'forecast_conditional_info',[],... 
        'forecast_regime',[],...
        'forecast_start_date','',...
        'forecast_end_date','',...
        'forecast_historical_db',rise_time_series.empty(0),...
        'truncated_normal_simulation_method','ghk');% alternative is 'gibbs'% should be moved elsewhere
    
    % OMGy,OMGx
    % alternatives for forecast_conditional_hypothesis are:
    % 'ncp'(or 'iris'): Number of conditioning periods
    % 'nas'(or 'iwb'): Number of anticipated steps
    % 0 (or 'jma'): No hypothesis
    
    return
end
%% pass the options
[obj.options,missing]=mysetfield(obj.options,varargin{:});
if ~isempty(missing)
    disp(missing)
    warning([mfilename,':: the options above are not recognized']) %#ok<WNTAG>
end
%% re-collect the options
forecast_param_uncert=obj.options.forecast_param_uncert;
forecast_shock_uncert=obj.options.forecast_shock_uncert;
forecast_state_uncert=obj.options.forecast_state_uncert;
if obj.markov_chains.regimes_number==1
    obj.options.forecast_regime=1;
end
forecast_regime_uncert=isempty(obj.options.forecast_regime);
forecast_conditional_info=obj.options.forecast_conditional_info;
forecast_historical_db=obj.options.forecast_historical_db;
forecast_start_date=obj.options.forecast_start_date;
forecast_end_date=obj.options.forecast_end_date;
forecast_paths=obj.options.forecast_paths;
% correct the number of paths in case there is no uncertainty
if ~forecast_shock_uncert && ~forecast_param_uncert
    disp([mfilename,':: no state uncertainty implemented yet'])
    forecast_paths=1;
end
truncated_normal_simulation_method=obj.options.truncated_normal_simulation_method;
forecast_conditional_hypothesis=obj.options.forecast_conditional_hypothesis;
%% Uncertainty in parameters
parameters=[];
if forecast_param_uncert
    posterior_density_ready=~isempty(vertcat(obj.estimated_parameters.mean));
    type='posterior_density';
    if ~posterior_density_ready
        mode_density_ready=~isempty([obj.estimation.mode]');
        type='mode_density';
        if ~mode_density_ready
            type='prior_density';
        end
    end
else
    parameters=vertcat(obj.estimated_parameters.mean);% 'mean'
    if isempty(parameters)
        parameters=[obj.estimation.mode]'; % 'mode'
        if isempty(parameters)
            parameters=[obj.estimation.priors.start]'; % 'calibration'
        end
    end
end

%% Uncertainty in shocks
%% Uncertainty in state or initial conditions
%% conditions

% if there are dates, then select the sample to pass to the smoother
% calibration change in smoothing and likelihood computation needs to be
% implemented

%% Load the data...
if isempty(forecast_historical_db)
    forecast_historical_db=obj.options.data;
    if isempty(forecast_historical_db)
        error([mfilename,':: historical data (initial conditions) not provided'])
    end
end
if forecast_historical_db.NumberOfPages~=1
    error('Historical database should have one page')
end
%% initial conditions
if isempty(forecast_start_date)
    smpl_end_date=forecast_historical_db.date_number(end);
    forecast_start_date=smpl_end_date+1;
else
    forecast_start_date=date2serial(forecast_start_date);
end
if isempty(forecast_end_date)
    smpl_end_date=forecast_historical_db.date_number(end);
    forecast_end_date=smpl_end_date+12;
else
    forecast_end_date=date2serial(forecast_end_date);
end

if forecast_end_date<=forecast_start_date
    error([mfilename,':: forecast_end_date<=forecast_start_date'])
end
histdb_varnames=forecast_historical_db.varnames;
model_varnames=obj.endogenous.name;
% are all the state variables present in the historical database? if not,
% the state is incomplete. This opens the way for state uncertainty
memb_locs=locate_variables(model_varnames,histdb_varnames,true);
incomplete_state=any(isnan(memb_locs));
forecast_state_uncert=forecast_state_uncert && incomplete_state;
% this is not enough, one should also check that observations are available
% for all variables at the origin of forecast.
start_vector=[];
if ~incomplete_state
    start_vector=forecast_historical_db(rise_date(forecast_start_date)-1);
    start_vector=start_vector(:);
    start_vector=start_vector(memb_locs);
    incomplete_state=any(isnan(start_vector));
    if incomplete_state
        start_vector=[];
    else
        start_vector=start_vector(:,ones(1,obj.markov_chains.regimes_number));
    end
end

if forecast_state_uncert||incomplete_state % we need to know where history starts
    if isempty(obj.options.estim_start_date)
        obj.options.estim_start_date=forecast_historical_db.start;
    end
    estim_start_date=obj.options.estim_start_date;
    estim_end_date=forecast_start_date-1;
    [obj,issue,retcode]=load_data(obj,'estim_start_date',estim_start_date,...
        'estim_end_date',estim_end_date,...
        'data',forecast_historical_db);
    disp(issue)
	if retcode
		error([mfilename,':: ',decipher_error(retcode)])
	end
end
%% conditional information

EndogenousConditions=[];
ShocksConditions=[];
if forecast_regime_uncert 
    if ~isempty(forecast_conditional_info)
        error([mfilename,':: conditional forecasting with multiple regimes not yet implemented. Try fixing a regime'])
    end
else
    regime_to_forecast=obj.options.forecast_regime;
    if ~ismember(obj.options.forecast_regime,(1:obj.markov_chains.regimes_number))
        error([mfilename,':: Chosen regime must be in [1,',int2str(obj.markov_chains.regimes_number),']'])
    end
    if ~isempty(forecast_conditional_info)
        endo_vars={};
        exo_vars={};
        fields=fieldnames(forecast_conditional_info);
        for ii=1:numel(fields)
            if ismember(fields{ii},{'lower_bound','upper_bound','central_tendency'})
                new_ts=forecast_conditional_info.(fields{ii});
                if new_ts.NumberOfPages~=1
                    error([mfilename,':: conditional database ',fields{ii},' must have one page only'])
                end
                if ii==1
                    min_date=new_ts.TimeInfo(1);
                    max_date=new_ts.TimeInfo(end);
                else
                    min_date=min(min_date,new_ts.TimeInfo(1));
                    max_date=max(max_date,new_ts.TimeInfo(end));
                end
                if ~isa(new_ts,'rise_time_series')
                    error([mfilename,'field ',fields{ii},' must be a rise_time_series object'])
                end
                if isempty(new_ts.varnames{1})
                    error([mfilename,':: conditional information variables should have names'])
                end
                for oo=1:new_ts.NumberOfVariables
                    vo=new_ts.varnames{oo};
                    if ismember(vo,obj.endogenous.name)
                        endo_vars=union(endo_vars,vo);
                    elseif ismember(vo,{obj.varexo.name}) 
                        exo_vars=union(exo_vars,vo);
                    else
                        error([mfilename,':: ',vo,' not recognized as an endogenous or an exogenous variable'])
                    end
                end
            else
                error([mfilename,':: conditional info fields must elements of lower_bound, upper_bound, central_tendency'])
            end
        end
        % check start date of the conditions with respect to the start of
        % the forecasts
        if min_date<forecast_start_date
            error([mfilename,':: conditional information cannot start before the beginning of the forecast period'])
        end
        addendum=numel(forecast_start_date:min_date)-1;
        date_span=min_date:max_date;
        date_numbers=[date_span.date_number];
        nobs=numel(date_span);
        all_vars=[endo_vars,exo_vars];
        restrictions_id=locate_variables(endo_vars,obj.endogenous.name);
        tmp=locate_variables(exo_vars,{obj.varexo.name});
        restrictions_id=[restrictions_id(:);tmp(:)];
        lower_bound=nan(nobs,numel(all_vars));
        upper_bound=nan(nobs,numel(all_vars));
        central_tendency=nan(nobs,numel(all_vars));
        for ii=1:numel(fields)
            new_ts=forecast_conditional_info.(fields{ii});
            dn=[new_ts.TimeInfo.date_number];
            date_locs=nan(size(dn));
            for oo=1:numel(dn)
                tmp=find(date_numbers==dn(oo));
                if isempty(tmp)
                    error([mfilename,':: date ',new_ts.TimeInfo(oo).date,' for ',fields{ii},' not found in the expanded date range'])
                end
                date_locs(oo)=tmp;
            end
            eval([fields{ii},'(date_locs,locate_variables(new_ts.varnames,all_vars))=double(new_ts);'])
        end
        % add the missing conditioning information dates...
        addendum=nan(addendum,numel(all_vars));
        lower_bound=[addendum;lower_bound];
        upper_bound=[addendum;upper_bound];
        central_tendency=[addendum;central_tendency];
        % holes are allowed in LB and UB only in places where the CT has
        % holes
        nan_ct=isnan(central_tendency);
        lower_bound(isnan(lower_bound) & ~nan_ct)=-inf;
        upper_bound(isnan(upper_bound) & ~nan_ct)=inf;
        % now format the conditions
        OMG_endo=[]; 
        OMG_exo=[]; % for the moment, I do not condition on a particular covariance matrix
        endo_locs=1:numel(endo_vars);
        exo_locs=endo_locs(end)+1:numel(restrictions_id);
        EndogenousConditions={central_tendency(:,endo_locs),OMG_endo,...
            struct('LB',lower_bound(:,endo_locs),'UB',upper_bound(:,endo_locs)),restrictions_id(endo_locs)};
        ShocksConditions={central_tendency(:,exo_locs),OMG_exo,...
            struct('LB',lower_bound(:,exo_locs),'UB',upper_bound(:,exo_locs)),restrictions_id(exo_locs)};
    end
end

%% do it

% number of forecasting steps
nsteps=numel(rise_date(forecast_start_date):rise_date(forecast_end_date));

history=1;
cond_fkst=nan(obj.endogenous.number(2),history+nsteps,forecast_paths);
cond_fkst_mean=nan(obj.endogenous.number(2),history+nsteps,forecast_paths);
uncond_fkst=nan(obj.endogenous.number(2),history+nsteps,forecast_paths);

ii=0;
PAI=1; % end-of-history probability in the constant case
while ii<forecast_paths
    % solve the model given the parameter draw
    if forecast_param_uncert
        % draw parameter
        parameters=obj.draw_parameter(type);
    else
        if ii==0 % solve the model only once
			obj=assign_estimates(obj,parameters);
			[obj,retcode]=solve(obj);
            if retcode
                error([mfilename,':: parameterization does not solve the model'])
            end
            parameters=[];
        end
    end
    
    [y0,retcode]=extract_initial_conditions();%obj,forecast_state_uncert
    
    if ~retcode
        ii=ii+1;
        % forecast_shock_uncert,forecast_regime,cond_data,cond_var_id,SS
        % extract the steady state from the SS_and_BGP
        SS=obj.steady_state_and_balanced_growth_path(1:obj.endogenous.number(2),:);
        if ii==1
            disp([mfilename,':: steady state problem unsolved'])
            disp([mfilename,':: switching steady states problem unsolved'])
            disp([mfilename,':: conditional forecast problem unsolved'])
            disp([mfilename,':: state uncertainty unsolved'])
            disp([mfilename,':: forecast_shock_uncert unsolved'])
            disp([mfilename,':: forecast_regime unsolved'])
        end
        if forecast_regime_uncert
            [cond_fkst(:,:,ii),cond_fkst_mean(:,:,ii),uncond_fkst(:,:,ii)]=multi_regime_forecast();
        else
            [cond_fkst(:,:,ii),cond_fkst_mean(:,:,ii),uncond_fkst(:,:,ii),junk,junk]=...
                forecast_engine(y0(:,regime_to_forecast),obj.T(:,:,regime_to_forecast),obj.R(:,:,:,regime_to_forecast),SS(:,regime_to_forecast));% compute the conditional forecasts
        end
    end
end

cond_fkst_db=rise_time_series(forecast_start_date-1,permute(cond_fkst,[2,1,3]),obj.endogenous.name);
cond_fkst_mean_db=rise_time_series(forecast_start_date-1,permute(cond_fkst_mean,[2,1,3]),obj.endogenous.name);
uncond_fkst_db=rise_time_series(forecast_start_date-1,permute(uncond_fkst,[2,1,3]),obj.endogenous.name);

% obj.forecasts=rise_time_series(obj.options.data.start,permute(y,[]),char(obj.endogenous.name));

% I need to select the conditional information. I still need to make up my
% mind on how to deal with the dichotomy of real_time vs soft conditions

% another issue

% and then splice the historical database and the forecast database
%
    function st=draw_state(probs)
        probs=cumsum(probs);
        probs(end)=1;
        st=find(probs>rand,1,'first');
    end

    function [mr_cond_fkst,mr_cond_fkst_mean,mr_uncond_fkst]=multi_regime_forecast()
        PAI_F=[PAI,nan(obj.markov_chains.regimes_number,nsteps)];
        mr_uncond_fkst=nan(obj.endogenous.number(2),nsteps+1,obj.markov_chains.regimes_number);
        mr_cond_fkst_mean=mr_uncond_fkst;
        for h=1:obj.markov_chains.regimes_number
            prev=h;
            mr_uncond_fkst(:,1,h)=y0(:,h);
            mr_cond_fkst_mean(:,1,h)=y0(:,h);
            for t=1:nsteps
                if h==1
                    PAI_F(:,t+1)=obj.Q'*PAI_F(:,t);
                end
                st=draw_state(PAI_F(:,t+1));
                mr_uncond_fkst(:,t+1,h)=SS(:,st)+obj.T(:,:,st)*(mr_uncond_fkst(:,t,h)-SS(:,prev));
                mr_cond_fkst_mean(:,t+1,h)=SS(:,st)+obj.T(:,:,st)*(mr_cond_fkst_mean(:,t,h)-SS(:,prev));
                if forecast_shock_uncert
                    mr_uncond_fkst(:,t+1,h)=mr_uncond_fkst(:,t+1,h)+obj.R(:,:,st)*randn(sum(obj.exogenous.number),1);
                end
                prev=st;
            end
        end
        mr_uncond_fkst=mean(mr_uncond_fkst,3);
        mr_cond_fkst_mean=mean(mr_cond_fkst_mean,3);
        mr_cond_fkst=mr_uncond_fkst;
    end

    function [y0,retcode]=extract_initial_conditions()%,forecast_historical_db,cond_db
        
        retcode=0;
        y0=start_vector;
        model_has_to_be_solved=~isempty(parameters);
        state_has_to_be_computed=isempty(y0);
        switch model_has_to_be_solved+2*state_has_to_be_computed
            case 0
            case 1
				obj=assign_estimates(obj,parameters);
                [obj,retcode]=solve(obj);
            case {2,3}
                [obj,junk,junk,retcode]=filter(obj,'evaluate_params',parameters);
                disp([mfilename,':: simulation smoother yet to be implemented'])
                if ~retcode
                    RegimesProbs=cell(1,numel(obj.markov_chains));
                    Regimes=obj.Regimes;
                    for kk=1:numel(obj.markov_chains)
                        reg_k=unique(Regimes(:,kk));
                        for istate=reg_k
                            data_k=double(obj.Filters.smoothed_probabilities.([obj.markov_chains(kk).name,'_',int2str(istate)]));
                            RegimesProbs{kk}=[RegimesProbs{kk},data_k(end)];
                        end
                    end
                    PAI=nan(obj.markov_chains.regimes_number,1);
                    for kk=1:obj.markov_chains.regimes_number
                        PAI(kk)=1;
                        for jj=1:numel(obj.markov_chains)
                            PAI(kk)=PAI(kk)*RegimesProbs{jj}(Regimes(kk,jj));
                        end
                    end
                    %
                    y0=nan(obj.endogenous.number(2),obj.markov_chains.regimes_number);
                    for kk=1:obj.endogenous.number(2)
                        vk=obj.varendo(kk).name;
                        data_k=double(obj.Filters.smoothed_variables.(vk));
                        y0(kk,:)=data_k(end,:);
                    end
                    clear data_k
                end
                if isempty(parameters)
                    % if the same parameter vector is to be used, then no
                    % need to recompute history, unless simulation
                    % smoothing is implemented. And this, I don't have time
                    % to do for now.
                    start_vector=y0;
                end
        end
        if obj.options.debug
            decipher_error(retcode)
        end
    end

    function [CYf,CYfMean,Yf,Emean,E]=forecast_engine(Y0,H,G,SS)
        % Description:
        % Usage:
        % Arguments:
        % Value:
        % Author(s): Junior Maih (junior.maih@norges-bank.no)
        % See Also:
        % Examples:
        NumberOfSimulations=1; % THIS IS TO REMAIN EXACTLY LIKE THIS !!!!
        
        %% demean the initial condition.
        yy0=Y0-SS;
        %% demean the endogenous conditions
        if ~isempty(EndogenousConditions)
            rest_id	=EndogenousConditions{4};
            SS_restr=transpose(SS(rest_id));
            EndogenousConditions{1}=bsxfun(@minus,EndogenousConditions{1},SS_restr);% [ncp,ncv]=size(MU);
            EndogenousConditions{3}.LB=bsxfun(@minus,EndogenousConditions{3}.LB,SS_restr);
            EndogenousConditions{3}.UB=bsxfun(@minus,EndogenousConditions{3}.UB,SS_restr);
        end
        %%
        [PbParams,M1,M2,RM2i,OMG,DYbar,junk,junk,Yf,MU,LB,UB]=...
            ConditionalProjectionSubEngine(forecast_conditional_hypothesis,nsteps,H,G,yy0,...
            EndogenousConditions,ShocksConditions);
        
        [rstar,rlow,rhigh,OMG]=ConditionalDistributionOfShocks(MU,OMG,LB,UB,DYbar);
        
        % check whether there are more holes
        if any(isnan(rlow))||any(isnan(rhigh))
            error([mfilename,':: LB and UB cannot have holes at other places than the holes in the Central tendency'])
        end
        if any(rlow>rhigh)
            error([mfilename,':: some elements in LB greater than their counterparts in UB'])
        end
        
        CSOMG=transpose(chol(OMG));
        % Compute the theoretical mean (doesn't change with the tightening in a symmetric
        % distribution)
        rmean=TruncatedMultivariateNormalMean(rstar,rlow,rhigh,CSOMG);
        gam2mean=RM2i*rmean;
        epsilon=M2*gam2mean;
        
        ShockHorizon=PbParams.NumberOfAnticipatedPeriods-1+PbParams.maxiter;
        Emean=InitializeArray(@zeros,PbParams.exo_nbr,ShockHorizon);
        Emean(:,1:PbParams.NumberOfCondShocksPeriods)=reshape(epsilon,PbParams.exo_nbr,[]);
        % Compute conditional forecasts
        CYfMean=InitializeArray(@zeros,PbParams.endo_nbr,PbParams.maxiter,1,yy0);
        CYfMean=ComputeForecasts(CYfMean,H,G,Emean,nsteps,PbParams.NumberOfAnticipatedPeriods);
        
        seed=[];
        r=TruncatedMultivariateNormalRnd(rstar,OMG,rlow,rhigh,...
            NumberOfSimulations,truncated_normal_simulation_method,seed);
        % draw gam2
        gam2=RM2i*r;
        clear r
        % draw gam1
        gam1=randn(PbParams.kk-PbParams.qq,NumberOfSimulations);
        % form epsilon
        epsilon=M1*gam1+M2*gam2;
        clear gam1 gam2
        
        % Build E
        E=InitializeArray(@randn,PbParams.exo_nbr,ShockHorizon,NumberOfSimulations);
        if ~isequal(forecast_conditional_hypothesis,0)
            E=0*E;
        end
        E(:,1:PbParams.NumberOfCondShocksPeriods,:)=reshape(epsilon,[PbParams.exo_nbr,PbParams.NumberOfCondShocksPeriods,NumberOfSimulations]);
        clear epsilon
        
        CYf=InitializeArray(@zeros,PbParams.endo_nbr,PbParams.maxiter,NumberOfSimulations,yy0);
        for isim=1:NumberOfSimulations
            % Compute conditional forecasts
            CYf(:,:,isim)=ComputeForecasts(CYf(:,:,isim),H,G,E(:,:,isim),nsteps,PbParams.NumberOfAnticipatedPeriods);
        end
        
        % Reset the forecasts to their normal length and re-mean them
        CYf=bsxfun(@plus,CYf(:,1:nsteps+1,:),SS);
        CYfMean=bsxfun(@plus,CYfMean(:,1:nsteps+1),SS);
        Yf=bsxfun(@plus,Yf(:,1:nsteps+1),SS);
        % and also the shocks (retain only the ones that matter)
        ShockHorizonFinal=PbParams.NumberOfAnticipatedPeriods-1+nsteps;
        E=E(:,1:ShockHorizonFinal,:);
        Emean=Emean(:,1:ShockHorizonFinal);
    end
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function [rstar,rlow,rhigh,OMG]=ConditionalDistributionOfShocks(MU,OMG,LB,UB,DTy)
rstar=MU-DTy;
rlow=LB-DTy;
rhigh=UB-DTy;

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function xx=TruncatedMultivariateNormalMean(m,lb,ub,CholCov)
npar=size(CholCov,1);
xx=nan(npar,1);
BOUNDS=[lb-m,ub-m];
for i=1:npar
    tmp=(BOUNDS(i,:)-CholCov(i,1:i-1)*xx(1:i-1,1))/CholCov(i,i);
    xx(i)=TruncatedNormalMean(tmp(1),tmp(2));
end
xx=m+CholCov*xx;
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function y=TruncatedNormalMean(lb,ub)
PHIl = NormalCumulativeDistribution(lb);
PHIr = NormalCumulativeDistribution(ub);
crit=1e-8;

finite_lb_and_ub=isfinite(lb) & isfinite(ub);
same=finite_lb_and_ub &  abs(ub-lb)<crit;
tails=abs(PHIr-PHIl)<crit & ~ same; 
good_tails=tails & finite_lb_and_ub;
bad_tails=tails & ~finite_lb_and_ub;
others=~same & ~tails; clear tails

if same % same, no problems
    y=ub;
elseif good_tails% assume a uniform distribution in the tails
    y=.5*(ub+lb);
elseif others% normal distribution for nice ones
    y=-(NormalDensity(ub)-NormalDensity(lb))/(PHIr-PHIl);
elseif bad_tails% Nasty ones
    y=0;
end
% if abs(lb-ub)<1e-9
%     y=ub;
% else
%     y=-(NormalDensity(ub)-NormalDensity(lb))/...
%         (NormalCumulativeDistribution(ub)-NormalCumulativeDistribution(lb));
% end
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function cdf=NormalCumulativeDistribution(x)
cdf=.5*(1+erf(x/sqrt(2)));
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function pdf=NormalDensity(x)
pdf=1/sqrt(2*pi)*exp(-.5*x^2);
end
