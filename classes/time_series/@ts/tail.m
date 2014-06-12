function db=tail(db,n)
if nargin<2
    n=min(5,db.NumberOfObservations);
end

span=fliplr(1:n)-1;
db=ts(db.date_numbers(end-span),db.data(end-span,:,:),db.varnames);

end