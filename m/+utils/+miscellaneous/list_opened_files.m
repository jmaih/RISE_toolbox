function X=list_opened_files()
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% lists all the files currently opened in the editor
X = matlab.desktop.editor.getAll;
X={X.Filename}';