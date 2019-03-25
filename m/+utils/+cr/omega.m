%  INTERNAL FUNCTION: creator of sparse matrices for mutivariate chain rules up to fifth order
% 
%  ::
% 
%    oo = omega(nz,index)
% 
%  Args:
% 
%     - **nz** [integer]: number of variables in the differentiation
%     - **index** [integer]: code for the requested omega. Must be in [1,9]
% 
%  Reference:
% 
%     - Oren Levintal(2014): "Fifth Order Perturbation Solution to
%       DSGE Models", 2014.
% 
%