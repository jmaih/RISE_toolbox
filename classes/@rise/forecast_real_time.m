function [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj,varargin)

% see also plot_real_time

default_options=struct('fkst_rt_nahead',16);

if isempty(obj)
    ts_fkst=default_options;
    return
end

nobj=numel(obj);
if nobj>1
    ts_fkst=cell(1,nobj);
    ts_rmse=cell(1,nobj);
    rmse=cell(1,nobj);
    Updates=cell(1,nobj);
    mystatenames={obj(1).varendo.name};
    for iobj=1:nobj
        if iobj>1 && ~isequal(mystatenames,{obj(iobj).varendo.name})
            error('models should have the same endogenous variables')
        end
        [ts_fkst{iobj},ts_rmse{iobj},rmse{iobj},Updates{iobj}]=...
            forecast_real_time(obj(iobj),varargin{:});
    end
    ts_rmse=concatenate_series_from_different_models(ts_rmse);
    % I probably should not do it this way since the unobservable variables
    % will never be the same.
    Updates=Updates{1}; 
    return
end

obj=set_options(obj,varargin{:});

fields=fieldnames(default_options);

for ifield=1:numel(fields)
    fi=fields{ifield};
    if isfield(obj.options,fi)
        default_options.(fi)=obj.options.(fi);
    end
end

nahead=default_options.fkst_rt_nahead;

[obj1,~,~,retcode]=filter(obj);

if retcode
    error(decipher_error(retcode))
end

Regimes=obj1.Regimes;
[nreg,nchains]=size(Regimes);
chain_names={obj1.markov_chains.name};
state_names=obj1.state_names;
state_probs=cell(1,numel(state_names));
for istate=1:numel(state_names)
    state_probs{istate}=double(obj1.Filters.updated_probabilities.(state_names{istate}));
end
regime_probs=cell(1,nreg);
for ireg=1:nreg
    regime_probs{ireg}=1;
    for ichain=1:nchains
        state=Regimes(ireg,ichain);
        loc=strcmp([chain_names{ichain},'_',int2str(state)],state_names);
        regime_probs{ireg}=regime_probs{ireg}.*state_probs{loc};
    end
end
regime_probs=transpose(cell2mat(regime_probs));
smpl=size(regime_probs,2);
endo_nbr=obj1.NumberOfEndogenous(2);
Updates=nan(smpl,endo_nbr,nreg);
state_var_names={obj.varendo.name};

for ivar=1:endo_nbr
    Updates(:,ivar,:)=double(obj1.Filters.updated_variables.(state_var_names{ivar}));
end
Updates=permute(Updates,[2,1,3]);

T=obj1.T;
Q=obj1.Q;

% R=obj1.R;

ss=vertcat(obj1.varendo.det_steady_state);

fkst=zeros(endo_nbr,nahead+1,smpl);

fkst_errors=nan(endo_nbr,nahead,smpl);

for t=1:smpl
    for ireg=1:nreg
        y0=Updates(:,t,ireg);
        yfkst=zeros(endo_nbr,nahead+1);
        yfkst(:,1)=y0;
        pai0=regime_probs(:,t);
        for h=1:nahead
            pai1=Q'*pai0;
            y1=iterate_forecast(y0);
            yfkst(:,1+h)=y1;
            y0=y1;
        end
         fkst(:,:,t)=fkst(:,:,t)+regime_probs(ireg,t)*yfkst;
    end
    hbar=min(nahead,smpl-t);
    steps=t+(1:hbar);
    fkst_errors(:,1:hbar,t)=Updates(:,steps)-fkst(:,1+(1:hbar),t);
end

rmse=nan(endo_nbr,nahead);

for h=1:nahead
    % discard the nans at the end
    this=squeeze(fkst_errors(:,h,1:end-h));
    rmse(:,h)=mean(this.^2,2);
end
rmse=sqrt(rmse);

ts_fkst=struct();
ts_rmse=struct();

% push the updates as a structure to ease plotting
Updates=obj1.Filters.Expected_updated_variables;

TimeInfo=Updates.(state_var_names{1}).TimeInfo;% [Updates.TimeInfo(2:end);Updates.TimeInfo(end)+1];
mynames=strcat(num2str((0:nahead)'),{'-step'});
for iendo=1:endo_nbr
    ts_fkst.(state_var_names{iendo})=...
        rise_time_series(TimeInfo,permute(fkst(iendo,:,:),[3,2,1]),mynames);
    ts_rmse.(state_var_names{iendo})=...
        rise_time_series(1,rmse(iendo,:)');
end

    function y1=iterate_forecast(y0)
        y1=0;
        for jreg=1:nreg
            y1j=one_step_forecast(jreg);
            y1=y1+pai1(jreg)*y1j;
        end
        function yout=one_step_forecast(index)
            yout=ss(:,index)+T(index)*(y0-ss(:,index));
        end
    end
end
