function data=bsxfun(db,fun,b)
% Overloaded bsxfun for ts object
%

data=bsxfun(fun,db.data,b);

data=ts(db.start,data);

end