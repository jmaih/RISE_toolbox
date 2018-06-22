function db=cumprod(db,dim)
% Overloaded cumprod function for ts object
%

if nargin<2
    dim=1;
end

db=ts(db.date_numbers,cumprod(db.data,dim),db.varnames);
end