%  myKnownRegimesFilter - filtering procedure for state-space models with
%  known switching dates
% 
%  ::
% 
% 
%    [LogLik,Incr,retcode,Filters]=myKnownRegimesFilter(syst,y,U,z,options,regs)
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
%     - **regs** [vector]: history of regimes, must be of the same length as
%        **y**
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
%     - If the filter is run on a constant parameter model, then the last
%       input argument need not be specified.
% 
%     - If the function is passed through a rise/dsge object, then it should
%       be called as
%       ff=filter(m,'kf_user_algo',{@myKnownRegimesFilter,regs}). That is,
%       RISE will provide all input arguments expect the last one that needs
%       to be provided by the user.
% 
%     - It is the responsibility of the user to make sure that the solution
%       of the model is consistent with agents expecting no transition from
%       one regime to another. This occurs when the original model is
%       backward-looking and/or the transition matrix is diagonal.
% 
%     See also: myConstantParamFilter
%