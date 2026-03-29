%--- help for rprt/rise_file ---
%
%  RISE_FILE Reads the content of a RISE model file and writes it to the report.
% 
%  Usage:
%    rise_file(obj, filename)
%    rise_file(obj, filename, options)
% 
%  Inputs:
%    - obj: The report object.
%    - filename: The name of the file to process.
%    - options (optional): Additional options for processing the file.
% 
%  Options:
%    - latexAlias: [true|false] Option to replace atoms with their LaTeX aliases. CURRENTLY NOT ACTIVE
%    - lineNumbers: [true|false] Add line numbers.
%    - lines: ['all'|empty|vector] Lines of code to report.
%    - syntax: [true|false] Highlight keywords
%    - footnote: [{''}|char] Add a footnote to the model.
% 
%  Example:
%    rise_file(obj, 'model_file.rs', 'lineNumbers', true);
%