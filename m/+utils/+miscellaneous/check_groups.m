%  check_groups : checks the adequacy of the partitions for any
%  decomposition
% 
%  Syntax
%  -------
%  ::
% 
%    [groups,positions]=check_groups(sList)
%    [groups,positions]=check_groups(sList,groups)
%  Inputs
%  -------
% 
%  - sList : [char|cell array] list of the variables(contributing factors);
% 
%  - groups : [structure|cell array |{empty}] grouping of shocks in the decomposition.
%    By default, the shocks are not grouped. The syntax is of the form
%    {group1,{v11,v12,...},...,groupn,{vn1,vn2,...}}. The shocks that are
%    not listed are put in a special group called "others". The "others"
%    group does not include the effect of initial conditions.
%    e.g. p=struct();
%         p.demand={'Ey','Er'};
%         p.supply={'Ep'};
%    e.g. p={'demand',{'Ey','Er'},'supply',{'Ep'}};
% 
%  Outputs
%  --------
% 
%  - groups : [struct] structure with the groupings. If the original groups
%    is empty, the decomposition is element by element. If there are
%    variables not included in the decomposition, a field called "others"
%    that gathers all of the unclassified variables.
% 
%  - positions : [cell array] each cell contains a vector indexes/locations
%    of the partition 
% 
%  Examples
%  ---------
% 
%  See also: dsge/historical_decomposition
%