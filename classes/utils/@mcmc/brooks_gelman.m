%--- help for mcmc/brooks_gelman ---
%
%  brooks_gelman computes the :cite:`BrooksGelman1998` convergence diagnostics, 
%  both the parameteric and the non-parameteric versions
% 
%  INPUTS:
% 
%    - this : mcmc object
% 
%    - do_plot : [true|{false}|empty]
% 
%  OUTPUTS:
% 
%  - myoutput  [struct] : Contains
%    
%        - UDIAG [by 6] double   
% 
%            - 1st column: length of total sequence interval
%            - 2nd column: sum of length of within sequence intervals; used to compute mean length of within sequence intervals
%            - 3nd column: within sequence variance
%            - 4nd column: sum of within sequence variances; used to compute mean within sequence variances
%            - 5nd column: within sequence kurtosis
%            - 6nd column: sum of within sequence kurtoses; used to compute mean within sequence kurtoses
% 
%        - MDIAG [by 6] double   
% 
%            - 1st column: length of total sequence interval
%            - 2nd column: sum of length of within sequence intervals; used to compute mean length of within sequence intervals
%            - 3nd column: within sequence variance
%            - 4nd column: sum of within sequence variances; used to compute mean within sequence variances
%            - 5nd column: within sequence kurtosis
%            - 6nd column: sum of within sequence kurtoses; used to compute mean within sequence kurtoses
% 
%        - h    [numeric|empty] : handle to the figures for the multivariate and univariate plots
%