function this3=plus(this1,this2)
if isa(this1,'rise_time_series')
    if isa(this2,'rise_time_series')
        [this1,this2]=intersect(this1,this2);
        if ~isequal(this1.NumberOfVariables,this2.NumberOfVariables)
            error([mfilename,':: datasets must have same number of columns'])
        end
        newdata=double(this1)+double(this2);
    elseif isa(this2,'double') && isscalar(this2)
        newdata=double(this1)+this2;
    else
        error([mfilename,':: plus operation undefined for this case'])
    end
    % Here it does not make sense to have names any more. But
    % all the same, perhaps I should have a function to rename
    % the series?
    this3=rise_time_series(this1.start,newdata);
elseif isa(this2,'rise_time_series') && isa(this1,'double') && isscalar(this1)
    newdata=this1+double(this2);
    % Here it does not make sense to have names any more. But
    % all the same, perhaps I should have a function to rename
    % the series?
    this3=rise_time_series(this2.start,newdata);
else
    error([mfilename,':: plus operation undefined for this case'])
end
end
