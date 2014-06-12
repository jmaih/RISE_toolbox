function [flag,year,period,freq,frequency]=is_serial(x)

[freq,frequency,year,period]=serial2frequency(x);

flag=~isempty(freq);

end