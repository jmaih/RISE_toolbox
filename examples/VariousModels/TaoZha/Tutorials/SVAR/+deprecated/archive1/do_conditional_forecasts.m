function [fkst,bands,hdl]=do_conditional_forecasts(m,db,conddb,draws,options)

% do_conditional_forecasts -- do conditional forecast for a SVAR model.
%
% ::
%
%
%   [fkst,bands]=do_conditional_forecasts(m,db,conddb)
%
%   [fkst,bands]=do_conditional_forecasts(m,db,conddb,draws)
%
%   [fkst,bands]=do_conditional_forecasts(m,db,conddb,draws,options)
%
% Args:
%
%    - **m** [svar]: svar model object (estimated or not)
%
%    - **db** [struct|ts]: time series (historical database)
%
%    - **conddb** [struct|ts]: time series (database with conditional information)
%
%    - **draws** [struct]: draws from the sampler (structure with fields "x"
%    and "f".
%
%    - **options** [empty|struct]: structure with possibly some of the
%    following fields
%      - **ndraws** [empty|{200}]: number of parameter draws
%      - **cbands** [empty|{[10,20,50,80,90]}]:
%      - **do_plot** [empty|true|{false}]: plot the conditional forecasts
%      - **param_uncertainty** [empty|{true}|false]: parameter uncertainty
%      - **shock_uncertainty** [empty|true|{false}]: shock uncertainty
%      - **nsteps** [empty|12]: number of forecast steps
%
% Returns:
%    :
%
%    - **fkst** [struct]: all forecasts for all parameter draws
%
%    - **bands** [struct]: confidence bands for all variables.
%
%    - **hdl** [empty|handle]: handle to the figure
%
% Note:
%
% Example:
%
%    See also:

if nargin<5
    options=[];
    if nargin<4
        draws=[];
    end
end
if isempty(options)
    options=struct();
end
defaults={
    'cbands',[10,20,50,80,90],@(x)isempty(x)||all(x>0 & x<100),'all elements in cbands must be in [0,100]'
    'do_plot',true,@(x)islogical(x),'do_plot should be true or false'
    'param_uncertainty',true,@(x)islogical(x),'param_uncertainty should be true or false'
    'shock_uncertainty',true,@(x)islogical(x),'shock_uncertainty should be true or false'
    'ndraws',200,@(x)isscalar(x) && isreal(x) && (floor(x)==ceil(x)) && x>0,'ndraws should be a positive integer'
    'nsteps',12,@(x)isscalar(x) && isreal(x) && (floor(x)==ceil(x)) && x>0,'nsteps should be a positive integer'
    };
options=parse_arguments(defaults,options);

cbands=options.cbands;
ndraws=options.ndraws;
do_plot=options.do_plot;
if isempty(draws)
    param_uncertainty=false;
else
    param_uncertainty=options.param_uncertainty;
end
shock_uncertainty=options.shock_uncertainty;

% Format conditioning information
%-----------------------------------
db=pages2struct(db);
conddb=pages2struct(conddb);
condvarnames=fieldnames(conddb);
for ivar=1:numel(condvarnames)
    vname=condvarnames{ivar};
    if isfield(db,vname)
        db.(vname)=cat(1,db.(vname),conddb.(vname));
    else
        db.(vname)=conddb.(vname);
    end
end
endo_vars=condvarnames(ismember(condvarnames,m.endogenous.name));
exo_vars=condvarnames(ismember(condvarnames,m.exogenous.name));
if ~isempty(condvarnames-endo_vars-exo_vars)
    disp(condvarnames-endo_vars-exo_vars)
    error('the variables above are not recognized as endogenous or exogenous')
end
% Now we fold the data for the conditioning variables : historical
% information on the first page 
%---------------------------------------------------------------------
history_end_date=m.options.estim_end_date;
alias=date2serial(history_end_date);
allvarnames=fieldnames(db);
for ivar=1:numel(allvarnames)
    vname=allvarnames{ivar};
    if ismember(vname,condvarnames)
        nail=find(db.(vname).date_numbers==alias);
        vals=double(db.(vname));
        histdata=vals(1:nail);
        conddata=vals(nail+1:end);
        nhist=numel(histdata);
        ncond=numel(conddata);
        vals=cat(3,histdata,nan(nhist,1,ncond));
        vals(end,1,2:end)=conddata;
        db.(vname)=ts(db.(vname).start,vals);
    end
end

% push the new dataset into the model, essentially replacing the old
% databse
%----------------------------------------------------------------------
m=set(m,'data',db);

% Also push some settings that will not change throughout the forecast
% exercise
%-----------------------------------------------------------------------
m=set(m,...
    'simul_history_end_date',history_end_date,...
    'forecast_cond_endo_vars',endo_vars,...
    'forecast_cond_exo_vars',exo_vars,...
    'forecast_nsteps',options.nsteps,...
    'forecast_to_time_series',false,...
    'forecast_shock_uncertainty',shock_uncertainty,...
    'simul_regime',1);

% loop over parameters and generate forecasts
%---------------------------------------------
for idraw=1:ndraws
    if param_uncertainty
        % push the draw inside the model, rather than getting out first
        [~,m]=draw_parameter(m,draws);
    end
    % compute a forecast
    mycast=forecast(m);
    if idraw==1
        start_date=mycast{2}{1,2};
        vnames=mycast{2}{2,2};
        nvars=numel(vnames);
        fkst=zeros(m.options.forecast_nsteps+1,nvars,ndraws);
    end
    fkst(:,:,idraw)=mycast{1};
end

% reconstruct time series as a sequence of pages
%------------------------------------------------
fkst=ts(start_date,fkst,vnames);

% compute the prctiles
%----------------------
bands=prctile(fkst,cbands);

% as a struct to separate out each variable
%-------------------------------------------
fkst=pages2struct(fkst);
bands=pages2struct(bands);

% add titles names to the bands but just for the model variables
%----------------------------------------------------------------
model_vars=union(m.endogenous.name,m.exogenous.name);
cbands_names=parser.create_state_list('ci',cbands);
bands.(model_vars{1})=ts(start_date,double(bands.(model_vars{1})),cbands_names);
prototype=bands.(model_vars{1});
for ivar=2:numel(model_vars)
    bands.(model_vars{ivar})=reset_data(prototype,...
        double(bands.(model_vars{ivar})));
end

hdl=[];

if do_plot
    figure('name','Conditional forecasts on ygap');
    for ivar=1:m.endogenous.number
        subplot(m.endogenous.number,1,ivar)
        vname=m.endogenous.name{ivar};
        plot(bands.(vname),'linewidth',2)
        title(vname)
        if ivar==3
            leg=legend(get(bands.(vname),'varnames'));
            set(leg,'interp','none','Orientation','horizontal',...
                'location','SouthOutside')
        end
    end
end
end
