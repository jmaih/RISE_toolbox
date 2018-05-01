function out=irf(obj,varargin)
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

if isempty(obj)
    
    mydefaults=irf@generic(obj);
    
    mydefaults=[mydefaults
        {'irf_anticipate',true,@(x)islogical(x),...
        'irf_anticipate must be true or false'}];
        
    if nargout
        
        out=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end

    
    return
    
end

obj=set(obj,varargin{:});

for iobj=1:numel(obj)
    
    obj(iobj).options.simul_anticipate_zero=~obj(iobj).options.irf_anticipate;
    
end

out=irf@generic(obj);

end