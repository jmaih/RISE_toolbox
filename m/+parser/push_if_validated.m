function dd=push_if_validated(val,testfunc,type,name_file_line)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<3
    name_file_line=[];
end
long_message=~isempty(name_file_line);
if isempty(testfunc)
    testfunc=@(x)~isnan(x);
end
if testfunc(val)
    dd=val;
else
    if long_message
        error(['wrong specification of ',type,' value for ',...
            name_file_line{1},' in ',name_file_line{2},' at line(s) ',...
            name_file_line{3}])
    else
        error(['wrong specification of ',type])
    end
end
end
