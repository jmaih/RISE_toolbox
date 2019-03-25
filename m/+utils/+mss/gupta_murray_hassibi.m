%  Checks for the Mean Square Stability using the Gupta, Murray and Hassibi algorithm
% 
%  ::
% 
%    C=gupta_murray_hassibi(T,Q)
%    C=gupta_murray_hassibi(T,Q,crit)
%    C=gupta_murray_hassibi(T,Q,crit,fast)
% 
%  Args:
% 
%     - **T** [1 x h cell array]: Autoregressive part of the solution in each
%       regime
%     - **Q** [h x h matrix]: Transition matrix
%     - **crit** [numeric]: tolerance criterion
%     - **fast** [false|{true}]: use the fastest algorithm
% 
%  Returns:
%     :
% 
%     - **flag** [true|false]: result of the investigation on whether there is
%       MSS or not.
%     - **max_eig** [numeric]: maximum eigenvalue
%     - **C** [matrix]: matrix in the criterion for checking stability
% 
%  See also: costa_fragoso_marques
% 
%  References:
%    - Vijay Gupta, Richard Murray and Babak Hassibi (2003):
%      "On the Control of Jump Linear Markov Systems with Markov State Estimation."
% 
%