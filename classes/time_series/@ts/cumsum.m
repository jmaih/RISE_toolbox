function db=cumsum(db,dim)
% INTERNAL FUNCTION
%

if nargin<2
    dim=1;
end

db=ts(db.date_numbers,cumsum(db.data,dim),db.varnames);
end