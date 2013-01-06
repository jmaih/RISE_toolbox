function [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj,varargin)

% see also plot_real_time

default_options=struct('fkst_rt_nahead',16);

if isempty(obj)
    ts_fkst=default_options;
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

Updates=obj1.Filters.Expected_updated_variables;
old_Updates=Updates;

T=obj1.T;

% R=obj1.R;

ss=vertcat(obj1.varendo.det_steady_state);

Updates=rise_time_series.collect(Updates);

state_var_names={obj.varendo.name};

locs=locate_variables(state_var_names,Updates.varnames);

dF=double(Updates);

endo_nbr=obj1.NumberOfEndogenous(2);

smpl=Updates.NumberOfObservations;

fkst=zeros(endo_nbr,nahead+1,smpl);

fkst_errors=nan(endo_nbr,nahead,smpl);

for t=1:smpl
    y0=transpose(dF(t,locs));
    fkst(:,1,t)=y0;
    for h=1:nahead
        y1=iterate_forecast(y0);
        fkst(:,1+h,t)=y1;
        y0=y1;
    end
    hbar=min(nahead,smpl-t);
    steps=t+(1:hbar);
    fkst_errors(:,1:hbar,t)=transpose(dF(steps,locs))-fkst(:,1+(1:hbar),t);
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
TimeInfo=Updates.TimeInfo;% [Updates.TimeInfo(2:end);Updates.TimeInfo(end)+1];
mynames=strcat(num2str((0:nahead)'),{'-step'});
for iendo=1:endo_nbr
    ts_fkst.(state_var_names{iendo})=...
        rise_time_series(TimeInfo,permute(fkst(iendo,:,:),[3,2,1]),mynames);
    ts_rmse.(state_var_names{iendo})=...
        rise_time_series(1,rmse(iendo,:)');
end

% push the updates as a structure to ease plotting
Updates=old_Updates;

    function y1=iterate_forecast(y0)
        y1=ss+T*(y0-ss);
    end
end
