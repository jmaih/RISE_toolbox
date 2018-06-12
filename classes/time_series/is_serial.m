function flag=is_serial(x)
% INTERNAL FUNCTION
%

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