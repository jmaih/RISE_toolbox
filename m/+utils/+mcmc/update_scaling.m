%  Updates the log of the scaling parameter for mcmc algorithms
% 
%  ::
% 
%    log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,n,xi3)
%    log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,n,xi3,c3)
%    log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,n,xi3,c3,c_range)
% 
%  Args:
% 
%     - **log_c** [numeric]: initial scaling
%     - **accept_ratio** [numeric]: acceptance rate
%     - **alpha_range** [interval|{[]}]: target acceptance range
%     - **fixed_scaling** [true|false]: Do not update the scale of the algorithm
%     - **n** [integer]: iteration
%     - **xi3** [numeric]: appears in gam3 = c3*(n+1)^(-xi3). must lie in
%       (0.5,1)
%     - **c3** [numeric|{1}]: appears in gam3 = c3*(n+1)^(-xi3);
%     - **c_range** [interval|{[sqrt(eps),100]}]: range of variation of the
%       scaling parameter (not its log!!!)
% 
%  Returns:
%     :
% 
%     - **log_c** [numeric]: updated scaling
% 
%  Note:
%     Adapts formulae 20 in Blazej Miasojedow, Eric Moulines and Matti
%     Vihola (2012): "Adaptive Parallel Tempering Algorithm"
% 
%