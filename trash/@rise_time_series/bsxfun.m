function data=bsxfun(db,fun,b)
if ~isvector(b)
    error('the third argument is expected to be a vector')
end
data=bsxfun(fun,double(db),b);
data=rise_time_series(db.start,data,db.varnames);
end