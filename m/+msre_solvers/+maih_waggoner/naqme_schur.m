%  naqme_schur : solves the system A*X^2+B*X+C=0
% 
%  ::
% 
%    [X]=naqme_schur(A,B,C,slvOpts)
%    [X,eigvals,retcode]=naqme_schur(...)
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
%      - **slvOpts.checkStab** [true|{false}] : check stability of the system
% 
%      - **slvOpts.xplosRoots** [{true}|false] : if all solutions are
%        computed we still can restrict ourselves to solutions that do not
%        involve explosive roots
% 
%      - **slvOpts.debug** [true|{false}] : slvOpts.debug or not
% 
%  Returns:
%     :
% 
%      - **X** [n x n x h x k array] : Solution set (k solutions)
% 
%      - **eigvals** [empty] : Eigenvalues (Not computed)
% 
%      - **retcode** [numeric] : 0 if there is no problem
% 
%      - **xtras** [struct] : information on the nature of the different
%        solutions
% 
%  See also :  udc
%