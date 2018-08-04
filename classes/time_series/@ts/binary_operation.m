function db=binary_operation(db1,db2,op_string)
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

operation=str2func(['@(x,y)bsxfun(@',op_string,',x,y)']);

if isa(db1,'ts')
    
    if isa(db2,'ts')
        
        [db1,db2]=intersect(db1,db2);
        %         if ~isequal(db1.NumberOfVariables,db2.NumberOfVariables)
        %         end
        try
            
            newdata=operation(db1.data,db2.data);
            
        catch
            
            error([mfilename,':: mismatch in sizes of datasets'])
            
        end
        
    elseif isa(db2,'double')
        
        if isscalar(db2)||isequal(size(db1.data),size(db2))
            
            newdata=operation(db1.data,db2);
            
        else
            
            error('wrong size of inputs')
            
        end
        
    else
        
        error([mfilename,':: plus operation undefined for this case'])
        
    end
    
    date_numbers=db1.date_numbers;
    
elseif isa(db2,'ts') && isa(db1,'double') && isscalar(db1)
    
    newdata=operation(db1,db2.data);
    
    date_numbers=db2.date_numbers;
    
else
    
    error([mfilename,':: plus operation undefined for this case'])
    
end

db=ts(date_numbers,newdata);

end
