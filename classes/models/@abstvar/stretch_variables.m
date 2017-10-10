function vbls=stretch_variables(obj,vbls)

if obj.is_panel
    
    nvars=numel(vbls);
    
    vnames=cell(1,nvars);
    
    for iii=1:nvars
        
        vnames{iii}=strcat(vbls{iii},'_',obj.members);
        
    end
    
    vbls=[vnames{:}];
    
end

end