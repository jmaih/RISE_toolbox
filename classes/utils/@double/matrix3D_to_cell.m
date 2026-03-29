%--- help for double/matrix3D_to_cell ---
%
%  matrix3D_to_cell Converts a 3D matrix array into a cell array of matrices.
% 
%  Args:
%    matrixArray : 3D matrix array where each slice represents a matrix
% 
%  Returns:
%    cellArray : Cell array where each cell contains a matrix from matrixArray
% 
%  Example:
%    matrixArray = rand(2, 3, 4);  % 4 matrices of size 2x3
%    cellArray = matrix3D_to_cell(matrixArray);
% 
%  See also: cell_2_3D_matrix
%
%    Documentation for double/matrix3D_to_cell
%       doc double/matrix3D_to_cell
%
%