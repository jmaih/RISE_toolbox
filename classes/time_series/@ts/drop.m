function this=drop(this,varargin)
% Drops the variable from the time series
%
% ::
%
%    db=drop(db,'var_name');
%
% Args:
%    db (ts object): time series object
%    varargin (string): names of the variables to drop
%
% Returns:
%    :
%
%    - db (ts object): time series object with corresponding variables dropped
%

survive=true(1,this.NumberOfVariables);
for ii=1:length(varargin)
    ids=locate_variables(varargin{ii},this.varnames);
    survive(ids)=false;
end
tmp_data=this.data;
this=ts(this.date_numbers,tmp_data(:,survive),this.varnames(survive));
end
