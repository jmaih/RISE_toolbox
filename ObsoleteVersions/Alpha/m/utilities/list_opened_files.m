function X=list_opened_files()
% lists all the files currently opened in the editor
X = matlab.desktop.editor.getAll;
X={X.Filename}';