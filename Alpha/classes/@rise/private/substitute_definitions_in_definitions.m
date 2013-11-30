function defcell=substitute_definitions_in_definitions(definitions)
% Get rid of definitions
defcell=cell(0,2);
for ii=1:definitions.number
    def_=definitions.shadow_dynamic{ii};
    equality=strfind(def_,'=');
    % get rid of the semicolon
    defcell=[defcell;{def_(1:equality-1),def_(equality+1:end-1)}]; %#ok<AGROW>
end
for ii=1:definitions.number
    for jj=ii+1:definitions.number
        defcell{jj,2}=strrep(defcell{jj,2},defcell{ii,1},['(',defcell{ii,2},')']);
    end
end
