function this=log(this,keep_name)
if nargin<2
    keep_name=false;
end
if keep_name
    this=rise_time_series(this.start,log(double(this)),this.varnames);
else
    this=rise_time_series(this.start,log(double(this)));
end
end
