%--- help for cell/max_discrepancy ---
%
%  `max_discrepancy` -- Calculate the maximum discrepancy between two cell arrays.
% 
%  This function calculates the maximum discrepancy between two cellures `cell1` and `cell2` by comparing the absolute differences of their corresponding fields. It returns the maximum discrepancy as a scalar value.
% 
%  Syntax:
%    o = max_discrepancy(cell1, cell2)
% 
%  INPUTS:
%  - `cell1` [cell]: The first cell array for comparison.
%  - `cell2` [cell]: The second cell array for comparison.
% 
%  OUTPUT:
%  - `o` [double]: The maximum discrepancy between the two cellures.
% 
%  EXAMPLE:
%  - Calculate the maximum discrepancy between two cellures `cellA` and `cellB`:
%    ```matlab
%    max_diff = max_discrepancy(cellA, cellB);
%    ```
% 
%  See also: cellfun.
%
%    Documentation for cell/max_discrepancy
%       doc cell/max_discrepancy
%
%    Other uses of max_discrepancy
%
%       dsge/max_discrepancy    struct/max_discrepancy
%