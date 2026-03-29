%--- help for dsge/parameters_to_file ---
%
%  PARAMETERS_TO_FILE Writes the parameters of a model object to a file.
% 
%    PARAMETERS_TO_FILE(model, fname) writes the parameters of the model object
%    (typically a DSGE model) to a file with the specified filename.
% 
%    Inputs:
%    - `model`: Model object containing the parameters.
%    - `fname`: File name to write the parameters to. The extension should be ".m".
% 
%    Example:
%        parameters_to_file(model, 'parameters_file.m')
% 
%    See also: get
%