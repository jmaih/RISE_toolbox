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

[obj1,~,~,retcode]=filter(obj,'kf_nsteps',nahead);

if retcode
    error(utils.error.decipher(retcode))
end

ts_fkst= obj1.filtering.rolling_multistep.Expected_filtered_variables;
ts_rmse= obj1.filtering.rolling_multistep.Expected_rmse;
rmse=double(ts_rmse);
Updates=obj1.filtering.Expected_updated_variables;


end
