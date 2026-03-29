%--- help for rprt/quote ---
%
%  QUOTE Adds a block quote to the report.
% 
%  Usage:
%    obj.quote('Text')
%    obj.quote('Text', 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Text: The text to include in the block quote.
%    - OptionName/OptionValue: Optional. Additional options for the block quote.
% 
%  Options:
%    - 'author': The author of the quote (default: '').
% 
%  Example:
%    rprt.quote('Lorem ipsum dolor sit amet.', 'author', 'John Doe');
%