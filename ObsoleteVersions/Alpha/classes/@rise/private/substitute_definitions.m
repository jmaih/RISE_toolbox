function string=substitute_definitions(string,defcell)
for ii=1:size(defcell,1)
    string=strrep(string,defcell{ii,1},['(',defcell{ii,2},')']);
end
end
