%--- help for rprt/subsection ---
%
%  SUBSECTION Adds a subsection to the report.
% 
%  Usage:
%    obj.subsection(Title)
%    obj.subsection(Title, 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Title: The title of the subsection.
%    - OptionName/OptionValue: Optional. Additional options for the subsection.
% 
%  Options:
%    - 'numbered': Boolean indicating whether the subsection should be numbered (default: true).
% 
%  Example:
%    rprt.subsection('Overview', 'numbered', false);
%    rprt.subsection('Methodology');
%