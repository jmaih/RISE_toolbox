%  INTERNAL FUNCTION: third-order multivariate chain rule
% 
%  ::
% 
%    res=third_order(dvvv,dvv,dv,vzzz,vzz,vz)
%    res=third_order(dvvv,dvv,dv,vzzz,vzz,vz,options)
% 
%  Args:
% 
%     - **dvvv** [nd x nv^3 matrix]: matrix of third derivatives of the d
%       function with respect to its locations. The derivatives are unfolded
%       columnwise
%     - **dvv** [nd x nv^2 matrix]: matrix of second derivatives of the d
%       function with respect to its locations. The derivatives are unfolded
%       columnwise
%     - **dv** [nd x nv matrix]: jacobian of function with respect to the
%       locations of its arguments
%     - **vzzz** [nv x nz^3 matrix]: third derivatives of the locations with
%       respect to the variables to differentiate. The derivatives are unfolded
%       columnwise
%     - **vzz** [nv x nz^2 matrix]: second derivatives (hessian) of the
%       locations with respect to the variables to differentiate. The derivatives
%       are unfolded columnwise
%     - **vz** [nv x nz matrix]: jacobian of the locations with respect to the
%       variables to differentiate
%     - **options** [empty|struct]: When not empty, options is a structure with
%       fields:
% 
%       - **large** [true|{false}] if true, a computation explicitly using the
%         kronecker product is avoided.
%       - **multiply** [true|{false}]: if true, explicit omega matrices are
%         constructed and then multiplied to other matrices to sum the
%         permutations. Else, a functional form is used instead.
% 
%  Returns:
%     :
% 
%     - **res** [nd x nz^3]: output matrix
% 
%