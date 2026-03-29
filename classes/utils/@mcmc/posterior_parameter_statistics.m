%--- help for mcmc/posterior_parameter_statistics ---
%
%  posterior_parameter_statistics computes posterior parameter statistics
% 
%  ::
% 
%     s = posterior_parameter_statistics(this)
%     s = posterior_parameter_statistics(this,Opts)
% 
%  Args:
% 
%     - this [mcmc object]: mcmc object
% 
%     - Opts [empty|struct]: options with the following fields
% 
%        - **percnt** [vector|{[2.5,25,50,75,97.5]}] : vector used in the
%          computation of percentiles
% 
%        -** mh_conf_sig** [scalar|{0.9}] : used in the computation of Highest
%          Probability Density Intervals.
% 
%        - **npoints** [integer|{512}] : number of points used  in the
%          computation of the density for each parameter
% 
%        - **kernel** [char|{'normal'}] : kernel used in the computation of the
%          density. Choices are 'epanechnikov'|'triangular'|'triweight',...
%          'uniform'|'cosine'
% 
%        - **do** [struct] : structure whose fields are the same as the ones for
%          the output. Each item should be true or false, with true triggering
%          the computation of the relevant item.
% 
%  Returns:
% 
%     **s** (struct): output with fields
% 
%        - posterior_mean [struct] : posterior mean
% 
%        - posterior_var [struct] : posterior variance
% 
%        - posterior_std [struct] : posterior standard deviation
% 
%        - posterior_percentiles [struct] : posterior posterior_percentiles
%          as defined by the user
% 
%        - posterior_deciles [struct] : posterior deciles 0.1:0.1:0.9
% 
%        - posterior_hpd_interval [struct] : Highest probability density intervals
% 
%        - posterior_density [struct] : kernel density estimate
%