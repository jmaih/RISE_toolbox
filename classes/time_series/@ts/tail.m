function db=tail(db,n)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    n=min(5,db.NumberOfObservations);
end

span=fliplr(1:n)-1;
db=ts(db.date_numbers(end-span),db.data(end-span,:,:),db.varnames);

end