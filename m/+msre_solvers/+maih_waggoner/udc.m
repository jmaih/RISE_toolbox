%  udc Undetermined coefficients solution algorithm for DSGE models.
%  The procedure can find all possible solutions for a constant-parameter
%  DSGE model or a regime-switching DSGE model with diagonal transition
%  matrix.
% 
%  ::
% 
%    [Tz_pb,eigvals,retcode]=udc(A,B,C,Q,T0,TolFun,maxiter)
% 
%  Args:
% 
%     - **A** [n x n] : Coefficient matrix on lead variables
% 
%      - **B** [n x n] : Coefficient matrix on current variables
% 
%      - **C** [n x n] : Coefficient matrix on lagged variables
% 
%      - **slvOpts.allSols** [true|{false}] : flag for finding all solutions
% 
%      - **slvOpts.msvOnly** [{true}|false] : return only MSV solutions
% 
%      - **slvOpts.debug** [true|{false}] : slvOpts.debug or not
% 
%  Returns:
%     :
% 
%      - **Tz_pb** [n x n x h x k array] : Solution set (k solutions)
% 
%      - **eigvals** [empty] : Eigenvalues (Not computed)
% 
%      - **retcode** [numeric] : 0 if there is no problem
% 
%  See also :  groebner
%