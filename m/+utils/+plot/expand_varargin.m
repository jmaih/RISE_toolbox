function vargs=expand_varargin(vargs,varargin)

if isempty(vargs)
    
    vargs={'LineWidth',2.5};
    
end

if isempty(varargin)
    
    varargin=vargs;
    
else
    
    found=false;
    
    for ii=1:2:length(varargin)-1
        
        found=strcmpi(varargin{ii},vargs{1});
        
        if found
            
            break
            
        end
        
    end
    
    if ~found
        
        varargin=[vargs,varargin];
        
    end
    
end

vargs=varargin;

end