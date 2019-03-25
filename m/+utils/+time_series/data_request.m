%  INTERNAL FUNCTION: Selects the data requested for estimation or for forecasting
% 
%  ::
% 
%     [data_values,start_date,end_date,missing]=data_request(data,varlist)
%     [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date)
%     [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date,end_date)
%     [data_values,start_date,end_date,missing]=data_request(data,varlist,start_date,end_date,pages)
% 
%  Args:
% 
%     - **data** [ts|structure]: database
%     - **varlist** [char|cellstr] : list of the variables of interest
%     - **start_date** [valid ts date|{''}]: start date of the sample requested
%     - **end_date** [valid ts date|{''}]: end date of the sample requested
%     - **pages** [integer|vector|{[]}]: pages requested
% 
%  Returns:
%     :
% 
%     - **data_values** [matrix|3-dimensional array]: data requested
%     - **start_date** [ts date]: start date of the sample requested
%     - **end_date** [ts date]: end date of the sample requested
%     - **missing** [integer]: number of missing observations replaced with
%        nans
% 
%  Note:
% 
%     - If there are insufficient data, the data are augmented with nans
% 
%