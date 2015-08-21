function names=create_names(names,prefix,n)
% CREATE_NAMES -- creates default variable names
%
% Syntax
% -------
% ::
%
%   names=CREATE_NAMES(names,prefix,n)
%
% Inputs
% -------
%
% - **names** [empty|char|cellstr]: if empty, names will be created. If
% not, they will be put to cellstr if necessary
%
% - **prefix** [char]: prefix of default names to be created. Must be a
% valid variable name
%
% - **n** [integer]: number of variable names to be created
%
% Outputs
% --------
%
% - **names** [cellstr]: created or transformed names
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if isempty(names)
    names=strcat({[prefix,'_']},cellstr(num2str((1:n)')));
    names=cellfun(@(x)x(~isspace(x)),names,'uniformOutput',false);
end
if ischar(names)
    names=cellstr(names);
end

end