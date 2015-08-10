function [pop,tags]=sort_population(pop)
% sort_population -- sorts population pf mcmc draws
%
% Syntax
% -------
% ::
%
%   [pop,tags]=sort_population(pop)
%
% Inputs
% -------
%
% - **pop** [struct]: vector of individuals, each with fields 
%   - **f** [numeric]: fitness level
%   - **x** [vector]: parameters 
%
% Outputs
% --------
%
% - **pop** [struct]: sorted vector of individuals
%
% - **tags** [vector]: reordering index
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

[~,tags]=sort([pop.f]);
pop=pop(tags);
end
