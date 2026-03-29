%--- help for rprt/equation ---
%
%  EQUATION Adds an equation to the report.
% 
%  Usage:
%    obj.equation('Expression')
%    obj.equation('Expression', 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Expression: The mathematical expression to display as an equation.
%    - OptionName/OptionValue: Optional. Additional options for the equation.
% 
%  Options:
%    - 'numbering': Boolean indicating whether the equation should be numbered (default: true).
% 
%  Example:
%    rprt.equation('E = mc^2', 'numbering', false);
%