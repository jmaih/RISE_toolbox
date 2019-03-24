%  Constructor for the time series (ts) object
% 
%  ::
% 
%    self=ts();    % construct a time series with no observations
%    self=ts(start_date,data);
%    self=ts(start_date,data,varnames);
%    self=ts(start_date,data,varnames,description);
%    self=ts(start_date,data,varnames,description,trailnans);
% 
%  Args:
% 
%     start_date (integer | char | serial date): start date of
%       the time series. The following are admitted:
% 
%       - annual data : e.g. 1990 or '1990'
%       - bi-annual data : e.g. '1990H1'
%       - Quarterly data : e.g. '1990Q3'
%       - monthly data : e.g. '1990M12'
% 
%     data (numeric): the format is nobs x nvars x npages,
%       where:
% 
%       - **nobs** is the number of observations
%       - **nvars** is the number of variables
%       - **npages** is the number of pages (3rd dimension)
% 
%     varnames (char | cellstr): names of the variables in the
%       database
%     description (char | cellstr | {''}): comments on each
%       variable in the database
%     trailnans (true|{false}): keep or remove nans (missing observations)
% 
%  Returns:
%     :
% 
%     - **self** [ts] : time series
% 
%
%    Reference page in Doc Center
%       doc ts
%
%