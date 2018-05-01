function [y,retcode]=exponential(x,g,c)
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
    
    c=0;
    
    if nargin<2
        
        error('at least 2 arguments should be provided')
        
    end
    
end

retcode=0;

y=nan;

if isa(g,'double') && any(g<0)
    
    if nargout>1
        
        retcode=4;
        
        return
        
    else
        
        error('g cannot be negative')
        
    end
    
end

y=1-exp(-g.*(x-c).^2);

end