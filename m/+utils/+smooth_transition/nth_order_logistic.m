function [y,retcode]=nth_order_logistic(x,g,varargin)
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
    
    error('at least 3 arguments should be provided')
    
end

n=length(varargin);

c1=varargin{1};

x_minus_c=x-c1;

retcode=0;

y=nan;

for ii=2:n
    
    c2=varargin{ii};
    
    if isa(c1,'double') && isa(c2,'double') &&  any(c2<c1)
        
        if nargout>1
            
            retcode=4;
            
            return
            
        else
            
            error('the c coefficients should be in an ascending order')
            
        end
        
    end
    
    x_minus_c=x_minus_c.*(x-c2);
    
    c1=c2;
    
end

if isa(g,'double') && any(g<0)
    
    if nargout>1
        
        retcode=4;
        
        return
        
    else
        
        error('g cannot be negative')
        
    end
    
end

y=1./(1+exp(-g.*x_minus_c));

end