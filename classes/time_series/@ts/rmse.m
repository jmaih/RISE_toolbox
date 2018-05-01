function [Rmse,Pe] = rmse(rawdata)
% RMSE -- Root Mean Square Error
%
% ::
%
%
%   [Rmse,Pe] = rmse(rawdata)
%
% Args:
%
%      - **rawdata** : [T x (h+1) x nsim ts]. Time series with
%          - T observations
%          - h+1 columns, where the first column represents the actual data
%          and the remaining h columns are forecasts
%          - nsim number of simulations
%
% Returns:
%    :
%
%      - **Rmse** : [h x nsim]. Matrix of root mean square errors
%
%      - **Pe** : [T x h x nsim ts]. Time series of prediction errors
%
% Note:
%
% Example:
%
%    See also:

% data = T x h x nsim

% the first column is the observed data

data=rawdata.data;

[T,hplus1,nsim]=size(data);

h=hplus1-1;

Pe=nan(T,h,nsim);

for t=1:T
    
    start=t+1;
    
    last=min(T,t+h);
    
    actual=data(start:last,1,:);
    
    if isempty(actual)
        
        continue
        
    end
    
    hstar=size(actual,1);
    
    predicted=data(t,1+(1:hstar),:);
    
    Pe(t,1:hstar,:)=permute(actual,[2,1,3])-predicted;
      
end

Mse = permute(nanmean(Pe.^2,1),[2,3,1]);

Rmse = sqrt(Mse);

Pe=ts(rawdata.start_date_number+1,Pe);

end




