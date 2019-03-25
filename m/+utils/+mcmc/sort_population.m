%  Sorts population pf mcmc draws
% 
%  ::
% 
%    [pop,tags]=sort_population(pop)
% 
%  Args:
% 
%     - **pop** [struct]: vector of individuals, each with fields
% 
%       - **f** [numeric]: fitness level
%       - **x** [vector]: parameters
% 
%  Returns:
%     :
% 
%     - **pop** [struct]: sorted vector of individuals
%     - **tags** [vector]: reordering index
% 
%