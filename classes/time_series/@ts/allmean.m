%--- help for ts/allmean ---
%
%  Compute all different type of means for the time series variables
% 
%  ::
% 
%     m = allmean(db);
% 
%  Args:
%     db (ts object): the times series object to compute the mean of
% 
%  Returns:
%     :
%     - m (cell): Cell containing
% 
%        - 1st column: strings denoting the type of mean (Harmonic, Geometric, Arithmetic, Quadratic)
%        - 1st row: variable names
%        - otherwise: mean value corresponding to the variable and mean concept.
% 
%