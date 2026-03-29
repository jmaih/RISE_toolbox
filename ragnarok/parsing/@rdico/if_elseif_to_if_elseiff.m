%--- help for rdico.if_elseif_to_if_elseiff ---
%
%  IF_ELSEIF_TO_IF_ELSEIFF - Replace all occurrences of "if_elseif" with "if_elseiff" in a function or function handle.
% 
%    out = if_elseif_to_if_elseiff(in)
% 
%  Input arguments:
%    - in: A function handle or a character vector representing a function.
% 
%  Output:
%    - out: The input function or function handle with "if_elseif" replaced by "if_elseiff".
% 
%  Description:
%  This function replaces all occurrences of "if_elseif" in the input function or function handle with "if_elseiff".
%  It can be used to update code that uses the deprecated "if_elseif" syntax to the modern "if_elseiff" syntax.
% 
%  Example:
%    % Define a function handle with "if_elseif" syntax
%    my_func = @(x) if_elseif(x < 0, -1, x > 0, 1, 0);
%    % Convert the function handle to use "if_elseiff" syntax
%    updated_func = if_elseif_to_if_elseiff(my_func);
% 
%    % Usage of the updated function handle
%    result = updated_func(5);
% 
%  Note: This function assumes that "if_elseif" and "if_elseiff" are used consistently and without nested function calls.
% 
%  See Also : utils.functions.replace_if_elseif
%