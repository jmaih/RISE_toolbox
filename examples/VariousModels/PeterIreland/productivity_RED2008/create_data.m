function data=create_data()

mydat=load('cih.dat');

data=ts('1948Q1',mydat,{'LG_C','LG_I','LG_H'});
end