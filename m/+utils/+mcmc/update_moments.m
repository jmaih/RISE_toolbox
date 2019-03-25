%  Updates the moments for mcmc algorithms
% 
%  ::
% 
%    [mu,SIG]=update_moments(mu,SIG,x,n,xi2)
%    [mu,SIG]=update_moments(mu,SIG,x,n,xi2,c2)
%    [mu,SIG]=update_moments(mu,SIG,x,n,xi2,c2,e)
% 
%  Args:
% 
%     - **mu** [vector]: mean
%     - **SIG** [matrix]: covariance
%     - **x** [vector|matrix]: new parameter vectors
%     - **n** [integer]: iteration
%     - **xi2** [numeric]: appears in gam2=c2*(n+1)^(-xi2). Must be in (0.5,1)
%     - **c2** [numeric|{1}]: appears in gam2=c2*(n+1)^(-xi2)
%     - **e** [numeric]: minimum eigenvalue of the updated covariance matrix
% 
%  Returns:
%     :
% 
%     - **mu** [vector]: updated mean
%     - **SIG** [matrix]: updated covariance matrix
% 
%  Note:
%     Applies formulae 23 and 24 in Blazej Miasojedow, Eric Moulines and Matti Vihola (2012): "Adaptive Parallel Tempering Algorithm"
% 
%