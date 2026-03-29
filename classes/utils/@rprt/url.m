%--- help for rprt/url ---
%
%  URL Adds a URL or hyperlink to the report.
% 
%  Usage:
%    obj.url('Text', 'Link')
%    obj.url('Text', 'Link', 'OptionName', OptionValue, ...)
% 
%  Inputs:
%    - obj: The report object.
%    - Text: The text of the URL.
%    - Link: The URL link.
%    - OptionName/OptionValue: Optional. Additional options for the URL.
% 
%  Options:
%    - 'newtab': Boolean indicating whether the URL should open in a new tab/window (default: false).
%    - MaskedText: Optional. The masked text to display instead of the URL (default: []).
% 
%  Example:
%    rprt.url('OpenAI', 'https://openai.com');
%    rprt.url('Click here', 'https://example.com', 'MaskedText','Website', 'newtab', true);
%