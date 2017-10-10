function self=set(self,varargin)

n=length(varargin);

if n
    
    tmp=reshape(varargin,2,[]);
    
    names=tmp(1,:);
    
    linres=strcmp(names,'linear_restrictions');
    
    varargin=reshape(tmp(:,~linres),1,[]);
    
    linres=tmp(:,linres);
    
    self=set@abstvar(varargin{:});
    
    if ~isempty(linres)
        
        if size(linres,2)~=1
            
            error('linear_restrictions appearing twice')
            
        end
        
        linres21=linres{2,1};
        
        if ischar(linres21)
            
            linres21=cellstr(linres21);
            
        end
        
        self.linear_restrictions_user=linres21;
        
        
    end
    
end

end
