%  INTERNAL FUNCTION: Efficient Kronecker Multiplication of kron(Q1,Q2,...,Qk)*X
% 
%  ::
% 
%    Y=kron_Q1_Qk_times_X(X,{Q1,n1},{Q2,n2},...,{Qk,nk})
% 
%  Args:
% 
%     - **X** [matrix|vector] is a vector or a matrix whose number of rows is equal to
%       the product of the columns of all Qi's
%     - **varargin** [array of 2-element cells] of the form {Qi,ni}
% 
%       - **Qi** [matrix]
%       - **ni** [integer] : number of times the kronecker of Qi appears
% 
%  Returns:
%     :
% 
%     - X
% 
%  Note:
% 
%     No kronecker product is required.
% 
%