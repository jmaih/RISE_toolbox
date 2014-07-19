function str=any2str(whatever)
if ischar(whatever)
    str=whatever;
elseif isnumeric(whatever)
    str=num2str(whatever);
elseif isa(whatever,'function_handle')
    str=func2str(whatever);
elseif iscell(whatever)
    str=parser.any2str(whatever{1});
else
    error([mfilename,':: cannot convert elements of the ',class(whatever),' class to a string'])
end
