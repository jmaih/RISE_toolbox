function db=plus(db1,db2)
% Overloaded + function for ts object
%

db=ts.binary_operation(db1,db2,mfilename);
end