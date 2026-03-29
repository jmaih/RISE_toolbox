%--- help for rprt/quotation ---
%
%  QUOTATION Adds a quotation block to the report.
% 
%  Usage:
%    obj.quotation('Text')
%    obj.quotation('Text', 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Text: The text to include in the quotation block.
%    - OptionName/OptionValue: Optional. Additional options for the quotation block.
% 
%  Options:
%    - 'author': The author of the quotation (default: '').
% 
%  Example:
%    rprt.quotation('Lorem ipsum dolor sit amet.', 'author', 'John Doe');
%