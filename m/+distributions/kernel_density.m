function [ff_,xx_]=kernel_density(data,lb,ub,kernel,n)
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

% kernel smoothing density estimation

kernel_functions={'epanechnikov',@(u)(3/4)*(1-u.^2).*(abs(u)<1)
    'normal',@(u)exp(-0.5 * u .^2) ./ (sqrt(2*pi))
    'triangular',@(u)(1-abs(u)).*(abs(u)<1)
    'triweight',@(u)(35/32)*((1-u.^2).^3).*(abs(u)<1)
    'uniform',@(u)(1/2)*(abs(u)<1)
    'cosine',@(u)pi/4*cos((pi/2)*u).*(abs(u)<1)
    };

if nargin<5
    
    n=[];
    
    if nargin<4
        
        kernel=[];
        
        if nargin<3
            
            ub=[];
            
            if nargin<2
                
                lb=[];
                
                if nargin<1
                    
                    if nargout>1
                        
                        error([mfilename,':: number of output arguments cannot exceed 1 when there are no inputs'])
                        
                    end
                    
                    disp(kernel_functions(:,1))
                    
                    return
                    
                end
                
            end
            
        end
        
    end
    
end

if isempty(kernel),kernel='normal';end

if isempty(n),n=100;end

if isempty(ub),ub=max(data);end

if isempty(lb),lb=min(data);end
    

if isempty(lb)
    
    lb=min(data);
    
end

if isempty(ub)
    
    ub=max(data);
    
end

% compute bandwith
h=std(data);

h=1.06*h*n^(-1/5); 

% normalizing function
normalize=@(xx,m,s)(xx-m)/s;

switch kernel
    
    case 'all'
        
        kk=1:size(kernel_functions,1);
        
    otherwise
        
        kk=find(strcmp(kernel,kernel_functions(:,1)));
        
end

if isempty(kk)
    
    error([mfilename,':: unrecognized kernel ',kernel])
    
end

xx=transpose(linspace(lb,ub,n));

ff=nan(n,numel(kk));

for ii=1:n
    
    u=normalize(xx(ii),data,h);
    
    for jj=1:numel(kk)
        
        Kfunc=kernel_functions{kk(jj),2};
        
        ff(ii,jj)=sum(Kfunc(u));
        
    end
    
end

ff=ff/(n*h);

if nargout==0
    
    plot(xx,ff)
    
    legend(kernel_functions(kk,1))
    
else
    
    xx_=xx;
    
    ff_=ff;
    
end



