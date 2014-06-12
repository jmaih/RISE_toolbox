function [year,period,freq,frequency]=date2year_period(date)
[serial_flag,year,period,freq,frequency]=is_serial(date);
if ~serial_flag
    date=date2serial(date);
    [serial_flag,year,period,freq,frequency]=is_serial(date);
    if ~serial_flag
        error('wrong date')
    end
end

end