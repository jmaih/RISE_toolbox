function data=bsxfun(db,fun,b)

data=bsxfun(fun,db.data,b);

data=ts(db.start,data);

end