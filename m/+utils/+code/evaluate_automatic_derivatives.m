%  INTERNAL FUNCTION: evaluate automatic differentiation
% 
%  ::
% 
%    C=evaluate_automatic_derivatives(funcs,order,tall_and_thin,engine,varargin)
% 
%  Args:
% 
%     - **funcs** [cell array]: functions to differentiate
%     - **order** [integer]: order of differentiation
%     - **engine** [@aplanar|@aplanar_]:
% 
%  Returns:
%     :
% 
%     - **C** [cell array]: derivatives
% 
%