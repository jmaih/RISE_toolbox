function db=times(db1,db2)
% Overloaded times function for ts object
%

db=ts.binary_operation(db1,db2,mfilename);
end