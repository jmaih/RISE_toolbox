%  returns a 3-column matrix in which:
%  - The first column represents the rows
%  - The second column represents the columns
%  - The third column represents the values expander, such that vals(col3)
%  returns a vector of same length as the rows of matrix B
% 
%  A is a cell array of cell arrays
%  - in each entry of A is a cell array
%    - where the first element is the row of the matrix
%    - the remaining entries are scalars or vectors that indicate the degree
%    of multiplicity of the corresponding derivative and the location of
%    those derivatives in the final matrix
%