function db=tail(db,n)
% Returns the last few sample dates of the time series
%
% ::
%
%    db = tail(db);
%    db = tail(db,n);
%
% Args:
%    db (ts object): time series object
%    n (integer): number of time steps to show
%
% Returns:
%    :
%    - db (ts object): time series with the last n-time steps
%
% Note:
%    - This is similar to the tail function in stata.
%

if nargin<2
    n=min(5,db.NumberOfObservations);
end

span=fliplr(1:n)-1;
db=ts(db.date_numbers(end-span),db.data(end-span,:,:),db.varnames);

end