function date_range=reset_range(date_range,dn)

if isempty(date_range)
    
    date_range=[dn(1),dn(end)];
    
end

if numel(date_range)~=2
    
    error('range must be a two-element vector')
    
end

if iscell(date_range)
    
    for iii=1:numel(date_range)
        
        if ischar(date_range{iii})
            
            date_range{iii}=date2serial(date_range{iii});
            
        end
        
    end
    
    date_range=cell2mat(date_range);
    
end

end