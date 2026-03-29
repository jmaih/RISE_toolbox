%  LOCATE_VARIABLES Locate variables in an array.
% 
%     IDs = LOCATE_VARIABLES(Variables, State, silent)
% 
%     Given a cell array of variables `Variables` and a cell array of state
%     variables `State`, this function returns the indices of the variables
%     in `Variables` in `State`.
% 
%     - `Variables`: Variable names to locate.
%     - `State`: Cell array of state variables.
%     - `silent` (optional): If true, suppress error messages for variables
%       not found. Default is false.
% 
%     Returns:
%     - `IDs`: Indices of the variables in `State`.
% 
%     Note:
%     - If `Variables` is empty, `IDs` will be an empty array.
%     - If `Variables` is a character array, it will be converted to a cell
%       array of strings.
%     - If `Variables` is a logical array, it should have the same number of
%       elements as `State`, and `IDs` will be the indices where the logical
%       array is true.
% 
%     Example:
%     - IDs = LOCATE_VARIABLES({'var1', 'var2'}, State);
% 
%     See also: OTHER_FUNCTION
%