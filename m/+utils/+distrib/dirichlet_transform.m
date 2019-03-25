%  INTERNAL FUNCTION: transforms the parameter of the dirichlet for estimation.
% 
%  ::
% 
%    aj=dirichlet_transform(x,sum_aij)
% 
%  Args:
% 
%     - **x** [k-1 x 1 vector]: vector of probabilities excluding the
%       "diagonal" element
%     - **sum_aij** [scalar]: sum of the weights including the "diagonal" element.
% 
%  Returns:
%     :
% 
%     - **aj** [k-1 x 1 vector]: un-normalized "weights" of the off-diagonal
%       elements
% 
%