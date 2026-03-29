%  Solves systems of non linear equations of several variables, using a
%  trust-region method.
% 
%  INPUTS
% 
%     - fcn:             objective function for which we need to find the
%       root
% 
%     - x0:              guess values
% 
%     - lb:              lower bound
% 
%     - ub:              upper bound
% 
%     - opts:            options (MaxIter,MaxFunEvals, etc.)
% 
%     - varargin:        further arguments of the fcn function
% 
%  OUTPUTS
% 
%     - x:               optimum
% 
%     - fval:            function value at the optimum
% 
%     - exitflag:        exit flag as in Matlab
% 
%     - output:          structure containing the details of the optimization
%       process
%