%--- help for ts/rmse ---
%
%  Compute the Root Mean Square Error
% 
%  ::
% 
%    [Rmse,Pe] = rmse(rawdata)
% 
%  Args:
% 
%       rawdata (T x (h+1) x nsim ts): Time series with
% 
%           - T observations
%           - h+1 columns, where the first column represents the actual data
%             and the remaining h columns are forecasts
%           - nsim number of simulations
% 
%  Returns:
%     :
% 
%       - **Rmse** : [h x nsim]. Matrix of root mean square errors
%       - **Pe** : [T x h x nsim ts]. Time series of prediction errors
% 
%