%--- help for rprt/subsubsection ---
%
%  SUBSUBSECTION Adds a subsubsection to the report.
% 
%  Usage:
%    obj.subsubsection(Title)
%    obj.subsubsection(Title, 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Title: The title of the subsubsection.
%    - OptionName/OptionValue: Optional. Additional options for the subsubsection.
% 
%  Options:
%    - 'numbered': Boolean indicating whether the subsubsection should be numbered (default: true).
% 
%  Example:
%    rprt.subsubsection('Overview', 'numbered', false);
%    rprt.subsubsection('Methodology');
%