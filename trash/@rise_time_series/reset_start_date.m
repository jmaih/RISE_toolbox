function this=reset_start_date(this,startdate)
if ~isequal(class(this),'rise_time_series')
    error([mfilename,':: input must be from class rise_time_series'])
end
this=rise_time_series(startdate,double(this),this.varnames);
end
