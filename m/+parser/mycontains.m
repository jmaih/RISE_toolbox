function flag=mycontains(str,substr)

try
    
    flag=contains(str,substr);
    
catch
    
    flag=~isempty(strfind(str,substr)); %#ok<STREMP>
    
end

end
