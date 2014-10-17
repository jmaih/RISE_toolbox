function [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj,varargin)
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


% see also plot_real_time

if isempty(obj)
    ts_fkst=struct('fkst_rt_nahead',16);
    return
end

nobj=numel(obj);
if nobj>1
    ts_fkst=cell(1,nobj);
    ts_rmse=cell(1,nobj);
    rmse=cell(1,nobj);
    Updates=cell(1,nobj);
    mystatenames=obj(1).endogenous.name;
    for iobj=1:nobj
        if iobj>1 && ~isequal(mystatenames,obj(iobj).endogenous.name)
            error('models should have the same endogenous variables')
        end
        [ts_fkst{iobj},ts_rmse{iobj},rmse{iobj},Updates{iobj}]=...
            forecast_real_time(obj(iobj),varargin{:});
    end
    ts_rmse=utils.time_series.concatenate_series_from_different_models(ts_rmse);
    % I probably should not do it this way since the unobservable variables
    % will never be the same.
    Updates=Updates{1}; 
    return
end

obj=set(obj,varargin{:});

nahead=obj.options.fkst_rt_nahead;

[obj1,~,~,retcode]=filter(obj);

if retcode
    error(utils.error.decipher(retcode))
end

nreg=obj1.markov_chains.regimes_number;

regime_probs=cell(1,nreg);
for ireg=1:nreg
    regime_probs{ireg}=double(obj1.filtering.smoothed_regime_probabilities.(['regime_',sprintf('%0.0f',ireg)]));
end
regime_probs=transpose(cell2mat(regime_probs));
smpl=size(regime_probs,2);
endo_nbr=obj1.endogenous.number(end);
Updates=nan(smpl,endo_nbr,nreg);
state_var_names=obj.endogenous.name;

for ivar=1:endo_nbr
    Updates(:,ivar,:)=double(obj1.filtering.updated_variables.(state_var_names{ivar}));
end
Updates=permute(Updates,[2,1,3]);

[T]=set_solution_to_companion(obj1);% T=obj1.solution.m_x;
Q=obj1.solution.transition_matrices.Q;% Q=obj1.solution.Q;

ss=obj1.solution.ss;

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
Updates=obj1.filtering.Expected_updated_variables;

TimeStart=Updates.(state_var_names{1}).start;% [Updates.TimeInfo(2:end);Updates.TimeInfo(end)+1];
mynames=strcat(num2str((0:nahead)'),{'-step'});
for iendo=1:endo_nbr
    ts_fkst.(state_var_names{iendo})=...
        ts(TimeStart,permute(fkst(iendo,:,:),[3,2,1]),mynames);
    ts_rmse.(state_var_names{iendo})=...
        ts(1,rmse(iendo,:)');
end

    function y1=iterate_forecast(y0)
        y1=0;
        for jreg=1:nreg
            y1j=one_step_forecast(jreg);
            y1=y1+pai1(jreg)*y1j;
        end
        function yout=one_step_forecast(index)
            yout=ss{index}+T{index}*(y0-ss{index});
        end
    end
end
