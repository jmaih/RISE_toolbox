function myfunc=remove_handles(myfunc)
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

if ischar(myfunc)
    
    myfunc=cellstr(myfunc);
    
end

if ~iscell(myfunc)
    
    error('input must be cell or char')
    
end

for ifunk=1:numel(myfunc)
    
    if isa(myfunc{ifunk},'function_handle')
        
        myfunc{ifunk}=func2str(myfunc{ifunk});
        
        first_closing=find(myfunc{ifunk}==')',1,'first');
        
        myfunc{ifunk}=myfunc{ifunk}(first_closing+1:end);
        
    end
    
end

end
