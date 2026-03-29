%--- help for mdd.global_options ---
%
%  GLOBAL_OPTIONS : generic options for the control of computations for
%  Marginal Data Density calculations. Key options include
% 
%    - **center_at_mean** [{false}|true]: if true, in the calculation of
%      moments, the covariance is centered at the mean, otherwise it is
%      centered at the mode. 
% 
%    - **L** [{500}|empty|integer]: number of extra draws
% 
%    - **debug** [{false}|true]: displays certain messages for debugging
%      purposes (algorithm specific)
% 
%    - **draws_are_iid** [{false}|true]:
% 
%    - and the options in fix_point_iterator()
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz
%