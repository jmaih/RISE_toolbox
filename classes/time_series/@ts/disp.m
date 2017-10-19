function disp(t,indent,fullChar)

if nargin < 3
    
    fullChar = false;
    
    if nargin < 2
        
        indent = 4;
        
    end
    
end

if isempty(t)
    
    disp('empty ts object')
    
    return
    
end

NumberOfVariables=t.NumberOfVariables;

NumberOfObservations=t.NumberOfObservations;

if ~((NumberOfObservations > 0) && (NumberOfVariables > 0))
    
    return
    
end

epilogue=index(t);

if NumberOfVariables==1
    
    epilogue=[t.description,{''},epilogue];
    
end

rownames=serial2date(t.date_numbers);

colnames=t.varnames;

the_data=t.data;

prologue=[];

table_displayer(the_data,colnames,rownames,prologue,epilogue,indent,fullChar)

end
