function [this1,this2]=intersect(this1,this2)
if nargin~=2
    error([mfilename,':: function must have 2 arguments'])
end
if ~isa(this1,'rise_time_series')||~isa(this2,'rise_time_series')
    error([mfilename,':: both arguments should be time series'])
end
if ~strcmp(this1.frequency,this2.frequency)
    error([mfilename,':: datasets must have same frequency'])
end
[~,I1,I2] = intersect(this1.date_number,this2.date_number);
if isempty(I1)||isempty(I2)
    error([mfilename,':: don''t have common dates'])
end
C=serial2date(this1.date_number(I1(1)));
vname1=this1.varnames;
vname2=this2.varnames;
this1=double(this1);
this2=double(this2);
this1=rise_time_series(C,this1(I1,:),vname1);
this2=rise_time_series(C,this2(I2,:),vname2);
end
