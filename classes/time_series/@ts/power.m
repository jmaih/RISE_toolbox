function db=power(db1,db2)
% Overloaded power function for ts object
%

db=ts.binary_operation(db1,db2,mfilename);
end