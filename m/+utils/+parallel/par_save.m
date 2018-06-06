function par_save(filename,variables,variables_names) %#ok<INUSL>
% INTERNAL FUNCTION
%

for ii=1:numel(variables_names)
    eval([variables_names{ii},'=variables{ii};'])
end
save(filename,variables_names{:})
