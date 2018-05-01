function X=list_opened_files()
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% lists all the files currently opened in the editor
X = matlab.desktop.editor.getAll;
X={X.Filename}';