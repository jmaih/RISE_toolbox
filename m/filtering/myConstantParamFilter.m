%  myConstantParamFilter - filtering procedure for state-space models with
%  constant parameters
% 
%  ::
% 
% 
%    [LogLik,Incr,retcode,Filters]=myConstantParamFilter(syst,y,U,z,options)
% 
%  Args:
% 
%     - **syst** [struct]: structure provided by dsge.filter
% 
%     - **y** [matrix]: matrix of data provided by dsge.filter
% 
%     - **U** [matrix]: matrix of trends provided by dsge.filter
% 
%     - **z** [matrix]: matrix of deterministic terms provided by dsge.filter
% 
%     - **options** [struct]: options provided by dsge.filter
% 
%  Returns:
%     :
% 
%     - **LogLik** [numeric]: value of the log likelihood
% 
%     - **Incr** [vector]: contributions to the likelihood in each period
% 
%     - **retcode** [numeric]: flag equal to 0 if there is no problem  
% 
%     - **Filters** [struct]: structure containing all the filtering
%       information
% 
%  Note:
% 
%     - If the function is passed through a rise/dsge object, then it should
%       be called as ff=filter(m,'kf_user_algo',@myConstantParamFilter).
% 
%     See also: myKnownRegimesFilter
%