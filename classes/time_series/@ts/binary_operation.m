function db=binary_operation(db1,db2,op_string)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

operation=str2func(['@(x,y)',op_string,'(x,y)']);
if isa(db1,'ts')
    if isa(db2,'ts')
        [db1,db2]=intersect(db1,db2);
        if ~isequal(db1.NumberOfVariables,db2.NumberOfVariables)
            error([mfilename,':: datasets must have same number of columns'])
        end
        newdata=operation(db1.data,db2.data);
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
