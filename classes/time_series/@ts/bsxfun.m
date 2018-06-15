function data=bsxfun(db,fun,b)
% INTERNAL FUNCTION
%

data=bsxfun(fun,db.data,b);

data=ts(db.start,data);

end