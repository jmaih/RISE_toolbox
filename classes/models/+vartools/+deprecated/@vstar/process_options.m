function [obj,varargin]=process_options(obj,varargin) 

n=length(varargin);

discard=false(1,n);

for iarg=1:2:n
    
    switch varargin{iarg}
            
%         case 'constant'
%             
%             obj.constant=varargin{iarg+1};
%             
%             if ~islogical(obj.constant)
%                 
%                 error('constant must be a logical')
%                 
%             end
%             
%             discard([iarg,iarg+1])=true;
%             
%         case 'nlags'
%             
%             obj.nlags=varargin{iarg+1};
%             
%             if ~(isreal(obj.nlags) && ...
%                     isnumeric(obj.nlags) && ...
%                     isscalar(obj.nlags) && ...
%                     obj.nlags>0 && ...
%                     floor(obj.nlags)==ceil(obj.nlags))
%                 
%                 error('nlags must be a positive constant')
%                 
%             end
%             
%             discard([iarg,iarg+1])=true;
            
    end
    
end

varargin(discard)=[];

end
