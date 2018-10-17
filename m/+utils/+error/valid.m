function flag=valid(x,check_real)
% INTERNAL FUNCTION: check validity of values of different data types
%

if nargin<2||isempty(check_real)
    
    check_real=true;
    
end

if iscell(x)

    flag=all(cellfun(@(x)utils.error.valid(x,check_real),x));

elseif isa(x,'double')
    
    flag=true;
    
    if check_real
        
        flag=isreal(x);
        
    end

    flag=flag && ~any(isnan(x(:))) && ~any(isinf(x(:)));

elseif isa(x,'tsparse')

    flag=utils.error.valid(x.v,check_real);

else

    error(['class ',class(x),' not a valid type for checking validity'])

end

end