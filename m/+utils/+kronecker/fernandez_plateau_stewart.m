%  INTERNAL FUNCTION: computes kron(Q1,Q2,...,Qk)*x or x*kron(Q1,...,Qk)
% 
%  ::
% 
%    Y=fernandez_plateau_stewart(do_left,X,{Q1,n1},{Q2,n2},...,{Qk,nk})
%    Y=fernandez_plateau_stewart(do_left,X,{Q1,n1},{Q2,n2},...,{Qk,nk},'-debug')
% 
%  Args:
% 
%     - **do_left** [true|false] if true, matrix x appears to the left of the
%       kronecker product i.e. x*kron(Q1,...,Qk). If instead do_left is false,
%       matrix x appears to the right of the kronecker product i.e.
%       kron(Q1,...,Qk)*x
%     - **x** [matrix|vector] is a vector or a matrix whose number of columns
%       (respectively number of rows) is equal to the product of the rows
%       (respectively columns) of all Qi's
%     - **varargin** [array of 2-element cells|matrix] cell inputs are of the
%       form {Qi,ni}. Alternatively, an input could simply be Qi
% 
%       - **Qi** [matrix] : must be square
%       - **ni** [integer] : number of times the kronecker of Qi appears
%       - **'-debug'** : entered at the end, this strings triggers the direct
%         computation of the x*kron(Q1,...,Qk) or kron(Q1,...,Qk)*x, checks the
%         differences and times each computation.
% 
%  Returns:
%     :
% 
%     - **x** [matrix|vector] : result of the computations
% 
%  Note:
% 
%     No kronecker product is required in the computation.
% 
%  Example:
% 
%     ::
% 
%        x=rand(1,3^3*7^2*5^1);
%        fernandez_plateau_stewart(true,x,{rand(3),3},{rand(7),2},{rand(5),1})
% 
%