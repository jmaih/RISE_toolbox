%--- help for cell/cell_2_3D_matrix ---
%
%  cell_2_3D_matrix Converts a cell array of matrices into a 3D matrix array.
% 
%  Args:
%    cellArray : Cell array where each cell contains a matrix of the same size
% 
%  Returns:
%    matrixArray : 3D matrix array where each slice represents a matrix from cellArray
% 
%  Example:
%    cellArray = {rand(2,3), rand(2,3), rand(2,3)};
%    matrixArray = cell_2_3D_matrix(cellArray);
% 
%  See also:
%
%    Documentation for cell/cell_2_3D_matrix
%       doc cell/cell_2_3D_matrix
%
%