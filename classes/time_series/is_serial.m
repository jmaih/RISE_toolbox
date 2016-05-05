function flag=is_serial(x)
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

if ~isnumeric(x)
    
    flag=false;
    
else
    
    d=x-floor(x);
    
    [~,unstamp]=time_frequency_stamp();
    
    freq=unstamp(d);
    
    flag=ismember(freq,[1,2,4,12,52]);
    
    flag=all(flag(:));
    
end

end