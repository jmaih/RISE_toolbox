% ./  Right array divide.
%    A./B divides each element of A by the corresponding element of B. A and 
%    B must have compatible sizes. In the simplest cases, they can be the 
%    same size or one can be a scalar. Two inputs have compatible sizes if, 
%    for every dimension, the dimension sizes of the inputs are either the 
%    same or one of them is 1. The element-wise operators ./ and .\ are 
%    related to each other by the equation A./B = B.\A.
% 
%    C = RDIVIDE(A,B) is called for the syntax 'A ./ B' when A or B is an object.
% 
%    See <a href="matlab:helpview('matlab','MATLAB_OPS')">MATLAB Operators and Special Characters</a> for more details.
% 
%    See also LDIVIDE, MLDIVIDE, MRDIVIDE.
%
%    Documentation for rdivide
%       doc rdivide
%
%    Other uses of rdivide
%
%       adolm/rdivide            rsymbdiff/rdivide    tabular/rdivide
%       codistributed/rdivide    splanar/rdivide      tall/rdivide
%       decomposition/rdivide    sym/rdivide          timeseries/rdivide
%       duration/rdivide         symbolic/rdivide     ts/rdivide
%       gpuArray/rdivide
%