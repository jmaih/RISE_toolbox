function [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj,varargin)
% forecast_real_time - forecast from each point in time
%
% ::
%
%
% - [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj)
% - [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj,varargin)
%
% Args:
%
%    - **obj** [dsge|svar|rfvar] : model object
%
%    - **varargin** : valid optional inputs coming in pairs. The main inputs
%      of interest for changing the default behavior are:
%      - **forecast_rt_nsteps** [integer] : number of periods ahead
%
% Returns:
%    :
%
%    - **ts_fkst** [struct] : fields are forecasts in the form of ts objects
%      for the different endogenous variables
%
%    - **ts_rmse** [ts|struct] : if only one object is processed, the output
%    is a ts. If instead several objects are processed, fields are RMSEs in
%    the form of ts objects for the different endogenous variables
%
%    - **rmse** [matrix] : RMSEs for the different endogenous variables
%
%    - **Updates** [struct] : fields are the updated (in a filtering sense) in
%      the form of ts objects for the different endogenous variables
%
% Note:
%
% Example:
%
%    See also: plot_real_time

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        ts_fkst=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
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
    
    ts_rmse=ts.concatenator(ts_rmse{:});
    
    % I probably should not do it this way since the unobservable variables
    % will never be the same.
    Updates=Updates{1}; 
    
    return
    
end

obj=set(obj,varargin{:});

nahead=obj.options.forecast_rt_nsteps;

[obj1,~,~,retcode]=filter(obj,'kf_nsteps',nahead);

if retcode
    
    error(utils.error.decipher(retcode))
    
end

ts_fkst= obj1.filtering.rolling_multistep.Expected_filtered_variables;

ts_rmse= obj1.filtering.rolling_multistep.Expected_rmse;

rmse=double(ts_rmse);

Updates=obj1.filtering.Expected_updated_variables;

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'forecast_rt_nsteps',16,@(x)num_fin_int(x) && x>=1,...
    'forecast_rt_nsteps must be a finite and positive integer'
    };

end
