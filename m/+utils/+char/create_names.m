%  INTERNAL FUNCTION: creates default variable names
% 
%  ::
% 
%    names=create_names(names,prefix,n)
% 
%  Args:
% 
%     - **names** [empty|char|cellstr]: if empty, names will be created. If
%       not, they will be put to cellstr if necessary
% 
%     - **prefix** [char]: prefix of default names to be created. Must be a
%       valid variable name
% 
%     - **n** [integer]: number of variable names to be created
% 
%  Returns:
%     :
% 
%     - **names** [cellstr]: created or transformed names
% 
%