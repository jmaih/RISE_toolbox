function this=drop(this,varargin)
survive=true(1,this.NumberOfVariables);
for ii=1:length(varargin)
    ids=locate_variables(varargin{ii},this.varnames);
    survive(ids)=false;
end
tmp_data=this.data;
this=ts(this.date_numbers,tmp_data(:,survive),this.varnames(survive));
end
