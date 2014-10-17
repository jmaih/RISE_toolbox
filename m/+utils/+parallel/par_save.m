function par_save(filename,variables,variables_names) %#ok<INUSL>
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

for ii=1:numel(variables_names)
    eval([variables_names{ii},'=variables{ii};'])
end
save(filename,variables_names{:})
