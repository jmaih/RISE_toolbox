function this=reset_start_date(this,startdate)

this=ts(startdate,this.data,this.varnames);
end
