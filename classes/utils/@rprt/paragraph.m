%--- help for rprt/paragraph ---
%
%  PARAGRAPH Adds a paragraph to the report.
% 
%  Usage:
%    obj.paragraph(Title)
%    obj.paragraph(Title, 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Title: The title of the paragraph.
%    - OptionName/OptionValue: Optional. Additional options for the paragraph.
% 
%  Options:
%    - 'numbering': Boolean indicating whether the paragraph should be numbered (default: true).
% 
%  Example:
%    rprt.paragraph('Overview');
%    rprt.paragraph('Methodology');
%