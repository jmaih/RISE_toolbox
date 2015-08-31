function [fkst,bands]=do_conditional_forecasts(m,db,draws,ndraws,cbands)
% do_conditional_forecasts -- do conditional forecast for a SVAR model.
%
% Syntax
% -------
% ::
%
%   [fkst,bands]=do_conditional_forecasts(m,db,draws)
%
%   [fkst,bands]=do_conditional_forecasts(m,db,draws,ndraws)
%
%   [fkst,bands]=do_conditional_forecasts(m,db,draws,ndraws,cbands)
%
% Inputs
% -------
%
% - **m** [svar]: svar model object (estimated or not)
%
% - **db** [struct|ts]: time series (historical database) 
%
% - **draws** [struct]: draws from the sampler (structure with fields "x"
% and "f".
%
% - **ndraws** [empty|{200}]: number of parameter draws 
%
% - **cbands** [empty|{[10,20,50,80,90]}]: 
%
% Outputs
% --------
%
% - **fkst** [struct]: all forecasts for all parameter draws
%
% - **bands** [struct]: confidence bands for all variables.
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if nargin<5
    cbands=[];
    if nargin<4
        ndraws=[];
    end
end
if isempty(cbands)
    cbands=[10,20,50,80,90];
end
if isempty(ndraws)
    ndraws=200;
end

% Generate conditions on ygap
%-----------------------------
nconds=14;
ygap=zeros(nconds,1);
ygap(1)=-0.025; % assumed value
for icond=2:nconds
    ygap(icond)=ygap(icond-1)+0.005;
end

% add those conditions as additional pages to the database
%---------------------------------------------------------
nobs=get(db.ygap,'NumberOfObservations');
tmp=double(db.ygap);
tmp=cat(3,tmp,nan(nobs,1,nconds));
tmp(end,1,2:end)=ygap;
db.ygap=ts(db.ygap.start,tmp);

% push the new dataset into the model, essentially replacing the old
% databse
%----------------------------------------------------------------------
m=set(m,'data',db);

% Also push some settings that will not change throughout the forecast
% exercise
%-----------------------------------------------------------------------
m=set(m,...
    'simul_history_end_date',db.ygap.finish,...
    'forecast_cond_endo_vars','ygap',...
    'forecast_nsteps',12,...
    'forecast_to_time_series',false,...
    'simul_regime',1);

% loop over parameters and generate forecasts
%---------------------------------------------
for idraw=1:ndraws
    % push the draw inside the model, rather than getting out first
    [~,m]=draw_parameter(m,draws);
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

end
