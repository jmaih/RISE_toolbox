%--- help for dsge/forecast_real_time ---
%
%  Forecast from each point in time
% 
%  ::
% 
%     [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj)
%     [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(obj,varargin)
% 
%  Args:
% 
%     obj (dsge | svar | rfvar): model object
%     varargin : valid optional inputs coming in pairs. The main inputs
%       of interest for changing the default behavior are:
% 
%       - **forecast_rt_nsteps** [integer] : number of periods ahead
% 
%  Returns:
%     :
% 
%     - **ts_fkst** [struct] : fields are forecasts in the form of ts objects
%       for the different endogenous variables
% 
%     - **ts_rmse** [ts|struct] : if only one object is processed, the output
%       is a ts. If instead several objects are processed, fields are RMSEs in
%       the form of ts objects for the different endogenous variables
% 
%     - **rmse** [matrix] : RMSEs for the different endogenous variables
% 
%     - **Updates** [struct] : fields are the updated (in a filtering sense) in
%       the form of ts objects for the different endogenous variables
% 
%  See also:
%     - plot_real_time
% 
%