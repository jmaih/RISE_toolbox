%--- help for dsge/predictive_analysis ---
%
%  `PREDICTIVE_ANALYSIS`: Perform predictive analysis for a dsge model
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(m, priors, userfun)
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,myDrawsMatrix)
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,'halton')
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,{'halton',N})
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,'sobol')
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,{'sobol',N})
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,'latin_hypercube')
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,{'latin_hypercube',N})
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,'prior')
% 
%    [draws, ok, retcode, timing, varargout] = predictive_analysis(...,{'prior',N})
% 
%  **Inputs**:
% 
%    - `m` (rise/dsge model object): RISE/DSGE model object.
%    - `priors` (structure): Structure with parameter names and their
%      distributions. 
%    - `userfun` (function_handle): Function that checks whether a
%      particular parameterization yields the desired outcome. 
%      The function should take as input the parameterized model (m) and
%      return at least two outputs in which the first is a boolean (true or
%      false) and the second is the retcode. The function can return
%      additional outputs that will be captured in varargout (see below).
% 
%    - `DrawsInfo`: could be
% 
%      - Matrix : this represents the draws computed/obtained elsewhere
%      - char : procedure to use for the draws ({'prior'}|'latin_hypercube'|'sobol'|'halton'):
%      - cell array : {procedure,N} where N is the number of draws with
%        default value 2^12.
% 
%  **Outputs**:
% 
%    - `draws` (numeric): Draws from the sampling.
%    - `ok` (logical): 1 x n vector of booleans, with true where a
%      particular parameter vector checks the desired behavior. 
%    - `retcode` (numeric): 1 x n vector of return codes, where 0 means
%      there are no problems. 
%    - `timing` (structure): Structure with the timing of various
%      sampling. 
%    - `varargout` (cell): User desired additional outputs beyond the two
%      first outputs of `userfun`
% 
%  **Notes**:
%    - The function performs prior predictive analysis by drawing parameter
%      values from specified distributions and checking their effects using a
%      user-defined function. 
%    - The function can exploit parallelization if workers are fired up. 
% 
%  **Example**:
% 
%    ```
%    [draws, ok, retcode, timing, o3,o4,...,on] = predictive_analysis(m, priors, userfun);
%    ```
% 
%  **See also**: `QUASI_MONTE_CARLO`
%