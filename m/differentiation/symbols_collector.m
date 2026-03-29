%  SYMBOLS_COLLECTOR - Collect symbolic variables used in a function or equation.
% 
%    [allsymb, input_list, ff, with_respect_to] = symbols_collector(ff, with_respect_to)
% 
%  Input arguments:
%    - ff: A function handle or cell array of function handles representing equations or functions.
%    - with_respect_to (optional): A cell array specifying variables with respect to which derivatives
%      will be taken.
% 
%  Output:
%    - allsymb: A cell array containing symbolic variables collected from the input functions or equations.
%    - input_list: A cell array representing input arguments or state variables.
%    - ff: A cell array of function handles (same as input) for equations or functions.
%    - with_respect_to: A cell array specifying variables with respect to which derivatives will be taken.
% 
%  Description:
%  This function collects symbolic variables used in the provided functions or equations. It extracts symbolic
%  variables from the equations or functions, combines them with input arguments or state variables (if provided),
%  and returns the complete list of symbolic variables. The optional with_respect_to argument can be used to specify
%  variables with respect to which derivatives will be taken, and these variables will be included in the output.
% 
%  Example:
%    % Define equations and variables
%    syms x y z;
%    eq1 = x + y^2;
%    eq2 = z - exp(x);
%    input_list = {'x', 'y', 'z'};
% 
%    % Collect symbolic variables
%    [allsymb, input_list, ff, with_respect_to] = symbols_collector({eq1, eq2}, {'x'});
% 
%    % The output variables contain the following:
%    % allsymb = {'x', 'y', 'z'}
%    % input_list = {'x', 'y', 'z'}
%    % ff = {eq1, eq2}
%    % with_respect_to = {'x'}
% 
%  Note: The function collects symbolic variables from the provided equations or functions and combines them
%  with input_list and with_respect_to to create the complete list of symbolic variables.
%