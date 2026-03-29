%--- help for rprt.code_formatter ---
%
%  CODE_FORMATTER - Returns a function handle that generates code for
%  formatting text in different formats.
% 
%    CODEFUNC = CODE_FORMATTER(ITEM, FRMT) returns a function handle
%    CODEFUNC that can be used to generate code for formatting text in
%    different formats. ITEM specifies the type of formatting, such as
%    'bold', 'italics', 'underline', or a color name like 'red', 'blue',
%    etc. FRMT specifies the target format, which can be 'pdf', 'latex', or
%    'html'.
% 
%    Example:
% 
%    codeFunc = code_formatter('bold', 'latex');
% 
%    latexCode = codeFunc('This text is bold'); % Generates LaTeX code:
%    '\textbf{This text is bold}'
% 
%    Supported formats:
% 
%    - 'pdf' or 'latex': Generates code in LaTeX format using LaTeX commands
%      (\textbf{}, \textit{}, \underline{}, \textcolor{}).
% 
%    - 'html': Generates code in HTML format using HTML tags (<b>, <i>, <u>, <span>).
% 
%    Note: The function supports a predefined set of colors for both LaTeX
%    and HTML formats. Additional colors can be added to the 'colors' cell
%    array as needed.
% 
%    See also: sprintf
%