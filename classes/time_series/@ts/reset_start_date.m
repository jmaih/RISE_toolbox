function this=reset_start_date(this,startdate)
% INTERNAL FUNCTION
%

this=ts(startdate,this.data,this.varnames);
end
