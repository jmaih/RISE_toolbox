%--- help for rsymbdiff.initialize_symbols ---
%
%  initialize_symbols static method
%  x = initialize_symbols(symvars, wrt) creates an array of rsymbdiff objects.
%  Inputs:
%    - symvars: Cell array of strings of the variables
%    - wrt: Subset of symvars listing the variables to differentiate
%  Outputs:
%    - x: Cell array of rsymbdiff objects representing the symbolic variables
%