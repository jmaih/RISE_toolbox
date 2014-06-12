function this=drop(this,varargin)
survive=true(1,this.NumberOfVariables);
for ii=1:length(varargin)
    ids=locate_variables(varargin{ii},this.varnames);
    survive(ids)=false;
end
tmp_data=double(this);
this=rise_time_series(this.start,tmp_data(:,survive),this.varnames(survive));
end
