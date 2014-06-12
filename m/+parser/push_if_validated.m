function dd=push_if_validated(val,testfunc,type,name_file_line)
if nargin<3
    name_file_line=[];
end
long_message=~isempty(name_file_line);
if ~testfunc(val)
    if long_message
        error(['wrong specification of ',type,' value for ',name_file_line{1},' in ',name_file_line{2},' at line(s) ',name_file_line{3}])
    else
        error(['wrong specification of ',type])
    end
end
dd=val;
end
