function [c]=get_reference(obj)% ,obj
if ~isempty(obj.ref)
    c=obj.ref;
elseif isnumeric(obj.func)
    c=sprintf('%0.16g',obj.func);
elseif isempty(obj.args)
    c=obj.func;
else
	disp([mfilename,':: this case should never happen, please contact junior.maih@gmail.com'])
    obj=commit(obj);
    c=obj.ref;
end

end