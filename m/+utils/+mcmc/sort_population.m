function [pop,tags]=sort_population(pop)
% sort_population -- sorts population pf mcmc draws
%
% ::
%
%
%   [pop,tags]=sort_population(pop)
%
% Args:
%
%    - **pop** [struct]: vector of individuals, each with fields
%      - **f** [numeric]: fitness level
%      - **x** [vector]: parameters
%
% Returns:
%    :
%
%    - **pop** [struct]: sorted vector of individuals
%
%    - **tags** [vector]: reordering index
%
% Note:
%
% Example:
%
%    See also:

[~,tags]=sort([pop.f]);

pop=pop(tags);

end
