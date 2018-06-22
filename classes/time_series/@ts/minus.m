function db=minus(db1,db2)
% Overloaded minus function for ts object
%

db=ts.binary_operation(db1,db2,mfilename);
end