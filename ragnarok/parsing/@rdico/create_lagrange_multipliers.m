%--- help for rdico.create_lagrange_multipliers ---
%
%  CREATE_LAGRANGE_MULTIPLIERS Creates Lagrange multipliers for optimization problems.
% 
%     mults = CREATE_LAGRANGE_MULTIPLIERS(n)
% 
%     mults = CREATE_LAGRANGE_MULTIPLIERS(n, instant)
% 
%     This function generates Lagrange multipliers for optimization problems. It
%     can create either a sequence of multipliers up to a certain index or just a
%     single multiplier.
% 
%     - `n`: Number of Lagrange multipliers to create or vector of indices
%       for which to create Lagrange multipliers
%     - `instant`: Logical flag indicating whether to create a single
%       multiplier or a sequence. Default is false.
% 
%     Returns:
%     - `mults`: Cell array of Lagrange multipliers.
% 
%     Example:
%     - mults = CREATE_LAGRANGE_MULTIPLIERS(5, false);
%     - mults = CREATE_LAGRANGE_MULTIPLIERS([28,29], true);
% 
%     See also:
%