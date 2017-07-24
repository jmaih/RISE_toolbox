function str=remove_description(str)

dblq=find(str=='"',2,'last');

while ~isempty(dblq)
    
    str(dblq(1):dblq(2))=[];
    
    dblq=find(str=='"',2,'last');
    
end

end