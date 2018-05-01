function [y,retcode]=second_order_logistic(x,g,c1,c2)
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

if nargin<4
    
    c2=[];
    
    if nargin<3
        
        c1=0;
        
        if nargin<2
            
            error('at least 2 arguments should be provided')
            
        end
        
    end
    
end

retcode=0;

y=nan;

if isempty(c2)
    
    c2=c1;
    
end

if isa(c1,'double') && isa(c2,'double') &&  any(c2<c1)
    
    if nargout>1
        
        retcode=4;
        
        return
        
    else
        
        error('c2 cannot be less than c1')
        
    end
    
end

if isa(g,'double') && any(g<0)
    
    if nargout>1
        
        retcode=4;
        
        return
        
    else
        
        error('g cannot be negative')
        
    end
    
end

y=1./(1+exp(-g.*(x-c1).*(x-c2)));

end