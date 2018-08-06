function db=head(db,n)
% Returns the first few sample dates of the time series
%
% ::
%
%    db = head(db);
%    db = head(db,n);
%
% Args:
%    db (ts object): time series object
%    n (integer): number of time steps to show
%
% Returns:
%    :
%    - db (ts object): time series with the first n-time steps
%
% Note:
%    - This is similar to the head function in stata.
%

if nargin<2
    
    n=min(5,db.NumberOfObservations);
    
end

db=ts(db.date_numbers(1:n),db.data(1:n,:,:),db.varnames);

end