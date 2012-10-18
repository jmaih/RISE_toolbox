function obj=set_properties(obj,varargin)
nargs=length(varargin);
if rem(nargs,2)~=0
    error([mfilename,':: arguments should enter by pairs'])
end
object_properties=properties(obj);
for ii=1:nargs/2
    propname=varargin{2*ii-1};
    if isempty(find(strcmp(propname,object_properties),1))
        error([mfilename,...
            ':: propname is not a valid property of class ',class(obj)])
    end
    propval=varargin{2*ii};
    obj.(propname)=propval;
end
end
