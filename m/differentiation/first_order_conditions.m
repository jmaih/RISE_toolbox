%  FIRST_ORDER_CONDITIONS Computes the first-order conditions of a system of equations.
% 
%  Syntax:
% 
%    DEqtns = first_order_conditions(symvars, EqtnsSym, wrtSymb, alien_list)
% 
%  Inputs:
% 
%    - symvars: A cell array of symbolic variables that appear in the equations.
%    - EqtnsSym: A symbolic expression or system of equations.
%    - wrtSymb: A symbolic variable or a cell array of symbolic variables with respect
%               to which the first-order conditions are computed.
%    - alien_list: (Optional) A list of functions unknown to the
%      differentiator
%    - debug: (Optional) inspect what the differentiator does
% 
%  Outputs:
%    - DEqtns: A cell array representing the first-order conditions with respect to
%              the specified symbolic variables.
% 
%  Description:
%    This function computes the first-order conditions of a system of equations with
%    respect to a set of symbolic variables. It uses symbolic differentiation to
%    compute the partial derivatives of the equations and constructs a cell array
%    containing the resulting first-order conditions.
% 
%  See also:
%    rsymbdiff.differentiate
%