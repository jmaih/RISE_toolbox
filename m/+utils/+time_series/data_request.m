function [data_values,start_date,end_date,missing]=data_request(...
data,varlist,start_date,end_date,pages)
% data_request - selects the data requested for estimation or for forecasting
%
% ::
%
%
%   - [data_values,start_date,end_date,missing]=data_request(data,varlist)
%   - [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date)
%   - [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date,end_date)
%   - [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date,end_date,pages)
%
% Args:
%
%    - **data** [ts|structure]: database
%    - **varlist** [char|cellstr] : list of the variables of interest
%    - **start_date** [valid ts date|{''}]: start date of the sample requested
%    - **end_date** [valid ts date|{''}]: end date of the sample requested
%    - **pages** [integer|vector|{[]}]: pages requested
%
% Returns:
%    :
%
%    - **data_values** [matrix|3-dimensional array]: data requested
%    - **start_date** [ts date]: start date of the sample requested
%    - **end_date** [ts date]: end date of the sample requested
%    - **missing** [integer]: number of missing observations replaced with
%       nans
%
% Note:
%
%    - If there are insufficient data, the data are augmented with nans
%
% Example:
%
%    See also:

if nargin<5
    
    pages=[];
    
    if nargin<4
        
        end_date='';
        
        if nargin<3
            
            start_date='';
        
        end
        
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

npages=The_data.NumberOfPages;

if isempty(pages)
    
    pages=1:npages;
    
end

if any(pages>npages)
    
    error('the anticipated horizon of shocks exceeds the number of "extra" pages in the dataset')

end

data_values=nan(nobs__,nvarobs,numel(pages));

for ivar=1:numel(ids)
    
    if isnan(ids(ivar))
        
        warning(['variable ',varlist{ivar},' not found in the database'])
        
        continue
        
    end
    
    data_values(:,ivar,pages)=DataValues(:,ids(ivar),pages);
    
end
% extend as necessary
%--------------------
sizDV=size(data_values);

data_values=[data_values;nan(missing,sizDV(2),numel(pages))];

data_values=permute(data_values,[2,1,3]);

% send forward only the data that are needed.
data_values=data_values(:,start:finish,:);

end