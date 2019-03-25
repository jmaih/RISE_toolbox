%  INTERNAL FUNCTION: sets the transformed dirichlet parameters back to probabilities
% 
%  ::
% 
%    x=dirichlet_untransform(aj,sum_aij)
% 
%  Args:
% 
%     - **aj** [k-1 x 1 vector]: un-normalized "weights" of the off-diagonal
%       elements
%     - **sum_aij** [scalar]: sum of "weights" including the "diagonal" element.
% 
%  Returns:
%     :
% 
%     - **x** [k-1 x 1 vector]: vector of probabilities excluding the
%       "diagonal" element
% 
%