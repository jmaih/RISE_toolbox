function db=cumprod(db,dim)
% INTERNAL FUNCTION
%

if nargin<2
    dim=1;
end

db=ts(db.date_numbers,cumprod(db.data,dim),db.varnames);
end