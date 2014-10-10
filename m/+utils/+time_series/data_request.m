function [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date,end_date)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<4
    end_date='';
    if nargin<3
        start_date='';
    end
end

The_data=ts.collect(data);
if isempty(start_date)
    start_date=The_data.start;
end
if isempty(end_date)
    end_date=The_data.finish;
end
start=date2obs(The_data.start,start_date);
finish=date2obs(The_data.start,end_date);
if start<=0
    error('estimation start date inferior to data start')
end
nobs__=The_data.NumberOfObservations;
missing=max(0,finish-nobs__);
if missing>0
    warning('estimation end date greater than data end date: nan observations will be added')
end
nvarobs=numel(varlist);
ids=locate_variables(varlist,The_data.varnames,true);
DataValues=double(The_data);
data_values=nan(nobs__,nvarobs);
for ivar=1:numel(ids)
    if isnan(ids(ivar))
        warning(['variable ',varlist{ivar},' not found in the database'])
        continue
    end
    data_values(:,ivar)=DataValues(:,ids(ivar));
end
% extend as necessary
%--------------------
sizDV=size(data_values);
data_values=[data_values;nan(missing,sizDV(2))];
data_values=transpose(data_values);
% send forward only the data that are needed.
data_values=data_values(:,start:finish);
