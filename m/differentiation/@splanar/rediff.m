%--- help for splanar.rediff ---
%
%  differentiate - differentiates vectors of splanar objects
% 
%  ::
% 
% 
%    derivs=splanar.differentiate(eqtns,nwrt)
%    derivs=splanar.differentiate(eqtns,nwrt,order)
%    derivs=splanar.differentiate(eqtns,nwrt,order,verbose)
% 
%  Args:
% 
%     - **eqtns** [splanar|cell array]: vector or cell array of splanar objects
% 
%     - **nwrt** [integer]: number of variables for which differentiation is
%       taken
% 
%     - **order** [integer|{1}]: order of differentiation
% 
%     - **alien_list** [empty|cellstr|char]: list of alien functions
%       (functions that RISE does not recognize and that are to be
%       differentiated by the user himself).
% 
%     - **verbose** [true|{false}]: displays information about the process e.g.
%       the amount of time it takes to differentiate each order
% 
%  Returns:
%     :
% 
%     - **derivs** [structure]: each element in the structure contains:
%       - **size** [2-colum vector]: number of rows and number of columns of
%         the compacted derivatives i.e. unique columns
%       - **derivatives** [vector of splanar]: derivatives
%       - **nwrt** [integer]: number of variables for which differentiation is
%         taken
%       - **order** [integer]: order of differentiation of the current
%         structure
%       - **nnz_derivs** [integer]: number of non-zero derivatives
%       - **partitions** [vector]: vector that help reconstruct the
%         expanded/grand derivative matrix.
%       - **map** []: empty field that will be used in a later stage
% 
%  Note:
% 
%  Example:
% 
%     See also:
%