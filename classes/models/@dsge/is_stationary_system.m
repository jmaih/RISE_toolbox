function flag=is_stationary_system(obj)
% IS_STATIONARY_SYSTEM -- checks whether a the model is stationary. i.e.
% does not containt trending variables
%
% ::
%
%
%   flag=IS_STATIONARY_SYSTEM(obj)
%
% Args:
%
%    - **obj**[rise|dsge]: model object
%
% Returns:
%    :
%
%    - **flag**[true|false]: true if the model is stationary
%
% Note:
%
%    - There is a difference between stability and stationarity
%      - stability refers to the system as a whole and conditions for
%      stability are often assessed through eigenvalues inside the unit circle
%      - stationarity refers to a scalar stochastic process. And such a
%      process will be said stationary if its first and second moment do not
%      vary with time.
%
% Example:
%
%    See also: RISE_GENERIC/IS_STABLE_SYSTEM

if isempty(obj)
    
    flag=cell(0,4);
    
    return
    
end

flag=isempty(get(obj,'endo_list(~stationary)'));

end