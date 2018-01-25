function mb=maxbytes()


if ispc
    
    ss=memory;
    
    mb=ss.MaxPossibleArrayBytes;
    
else
    
    [~,m]=unix('vm_stat | grep free');
    
    dblq=find(m==':');
    
    m=m(dblq+1:end);
    
    mb=str2double(m(~isspace(m)));
    
    mb=mb*4096;
        
end

end