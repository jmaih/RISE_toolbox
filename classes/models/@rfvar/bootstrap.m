%--- help for rfvar/bootstrap ---
%
%  Bootstrap the parameter values resampling :math:`\varepsilon_t` terms from the data sample
% 
%  ::
% 
%     RepsRun = bootstrap(rfvar, boot)
% 
%  Args:
%     rfvar (rfvar object): rfvar object
%     boot (integer): number of bootstrap simulations (default: 1000)
% 
%  Returns:
%     :
% 
%     - **RepsRun** : bootstrapped values of the estimated parameters
% 
%        - dimension 1: parameter levels
%        - dimension 2: different runs of bootstrap samples
% 
%  Example:
%     For example, one can compute variance of the estimate of the parameters by computing::
% 
%        RepsRun = bootstrap(rfvar, 10000);
%        var(RepsRun, [], 2)
% 
%