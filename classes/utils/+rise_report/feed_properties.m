function obj=feed_properties(classname,obj,varargin)
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


n=length(varargin);
if rem(n,2)
    error('arguments must come in pairs')
end
if n
    allprops=properties(rise_report.(classname));
    props=varargin(1:2:n-1);
    bad=~ismember(props,allprops);
    if any(bad)
        disp(props(bad))
        error(['the properties above are not valid properties for class:: ',mfilename])
    end
    vals=varargin(2:2:n);
    for ii=1:0.5*n
        obj.(props{ii})=vals{ii};
    end
end

end