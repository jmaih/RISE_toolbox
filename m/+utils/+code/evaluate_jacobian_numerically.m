%  INTERNAL FUNCTION: numerical evaluation of the jacobian of the objective function
% 
%  ::
% 
%    J=evaluate_jacobian_numerically(funcs,y,x,ss,param,def,s0,s1)
% 
%  Args:
% 
%     - **funcs** [fhandle|cell array]: function or functions to be
%       differentiated
%     - **y** [vector]: values of endogenous variables
%     - **x** [vector]: values of exogenous variables
%     - **ss** [vector]: steady state
%     - **param** [vector]: parameter vector
%     - **def** [vector]: values of definitions
%     - **s0** [scalar]: state today
%     - **s1** [scalar]: state tomorrow
% 
%  Returns:
%     :
% 
%     - **J** [matrix]: Numerical Jacobian of **funcs** at [y,x,sparam]
% 
%  Note:
%     It is assumed that the inputs are y,x,ss,param,sparam,def,s0,s1 as
%     ordered in parser.input_list()
% 
%