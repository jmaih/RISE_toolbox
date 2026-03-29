%--- help for estimable/print_estimation_results ---
%
%  `PRINT_ESTIMATION_RESULTS`: Display the results of estimation.
% 
%    retcode = print_estimation_results(obj)
%    retcode = print_estimation_results(obj, detail)
% 
%  **Args**:
%    - `obj` (estimable): Model object.
%    - `detail` (true|{false}): If true, the description of parameters is
%      given alongside the code name. If they are the same, i.e., the
%      description has not been provided, then the code name is given only once. 
% 
%  **Returns**:
%    - `retcode` (numeric): 0 if there is no problem.
% 
%  **Notes**:
% 
%  - If there are multiple objects in the array, each object's results will
%    be displayed. 
%  - The function displays information such as log-posterior,
%    log-likelihood, log-prior, etc. 
%  - The displayed results include information about the estimation sample,
%    solution algorithm, estimation algorithm, number of estimated parameters, 
%    number of function evaluations, and the time taken for estimation.
%  - If there are any issues, a list of issues is displayed.
% 
%  **Example**:
%    print_estimation_results(obj)
%    print_estimation_results(obj, true)
% 
%  **See also**: `TABLE_DISPLAYER`
%