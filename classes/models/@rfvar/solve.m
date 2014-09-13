function [obj,retcode]=solve(obj,varargin)%

if isempty(obj)
    obj=struct();
    return
end

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end
nobj=numel(obj);
if nobj>1
    retcode=nan(1,nobj);
    for iobj=1:nobj
        [obj(iobj),retcode(iobj)]=solve(obj(iobj));
    end
    return
end

[obj,retcode]=solve@svar(obj,varargin);

if ~retcode        
    obj=structural_form(obj);
end

