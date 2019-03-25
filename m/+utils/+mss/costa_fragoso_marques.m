%  checks for Mean Square Stability using the Costa, Fragoso and Marques algorithm.
% 
%  ::
% 
%    C=costa_fragoso_marques(T,Q)
%    C=costa_fragoso_marques(T,Q,crit)
% 
%  Args:
% 
%     - **T** [1 x h cell array]: Autoregressive part of the solution in each
%       regime
%     - **Q** [h x h matrix]: Transition matrix
%     - **crit** [numeric]: tolerance criterion
% 
%  Returns:
%     :
% 
%     - **flag** [true|false]: result of the investigation on whether there is
%       MSS or not.
%     - **max_eig** [numeric]: maximum eigenvalue
%     - **C** [matrix]: matrix in the criterion for checking stability
% 
%  See also: gupta_murray_hassibi
% 
%  References:
% 
%    - O.L.V. Costa, M.D. Fragoso and R.P. Marques (2004):
%      "Discrete-Time Markov Jump Linear Systems"
% 
%