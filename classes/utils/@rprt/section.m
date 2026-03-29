%--- help for rprt/section ---
%
%  SECTION Adds a section to the report.
% 
%  Usage:
%    obj.section(Title)
%    obj.section(Title, 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Title: The title of the section.
%    - OptionName/OptionValue: Optional. Additional options for the section.
% 
%  Options:
%    - 'numbered': Boolean indicating whether the section should be numbered (default: true).
% 
%  Example:
%    rprt.section('Introduction', 'numbered', false);
%    rprt.section('Data Analysis');
%