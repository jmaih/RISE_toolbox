function db=rdivide(db1,db2)
% Overloaded rdivide function for ts object
%

db=ts.binary_operation(db1,db2,mfilename);
end