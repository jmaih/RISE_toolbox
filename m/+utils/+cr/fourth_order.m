%  INTERNAL FUNCTION: fourth-order multivariate chain rule
% 
%  ::
% 
%    res=fourth_order(dvvvv,dvvv,dvv,dv,vzzzz,vzzz,vzz,vz)
%    res=fourth_order(dvvvv,dvvv,dvv,dv,vzzzz,vzzz,vzz,vz,options)
% 
%  Args:
% 
%     - **dvvvv** [nd x nv^4 matrix]: matrix of fourth derivatives of the d
%       function with respect to its locations. The derivatives are unfolded
%       columnwise
%     - **dvvv** [nd x nv^3 matrix]: matrix of third derivatives of the d
%       function with respect to its locations. The derivatives are unfolded
%       columnwise
%     - **dvv** [nd x nv^2 matrix]: matrix of second derivatives of the d
%       function with respect to its locations. The derivatives are unfolded
%       columnwise
%     - **dv** [nd x nv matrix]: jacobian of function with respect to the
%       locations of its arguments
% 
%     - **vzzzz** [nv x nz^4 matrix]: fourth derivatives of the locations with
%       respect to the variables to differentiate. The derivatives are unfolded
%       columnwise
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
%     - **res** [nd x nz^4]: output matrix
% 
%