function obj=swap_references(obj,new_ref,order)
old_ref=obj.ref;
if ~isempty(old_ref) && old_ref(1)~='G'
% if ~isempty(old_ref) && ~strncmp(old_ref,'G',1)
    underscores=strfind(old_ref,'_');
    if order==str2double(old_ref(underscores(1)+1:underscores(2)-1))
        obj.ref=new_ref;
    end
end
end
