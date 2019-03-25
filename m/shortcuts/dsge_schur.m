%  dsge_schur : Schur solution algorithm for DSGE models.
%  The procedure can find all possible solutions for a constant-parameter
%  DSGE model or a regime-switching DSGE model with diagonal transition
%  matrix.
% 
%  ::
% 
%    [Tz_pb,eigvals,retcode]=dsge_schur(Alead,Acurr,Alag,Q,T0,TolFun,maxiter)
% 
%    [Tz_pb,eigvals,retcode]=dsge_schur(Alead,Acurr,Alag,Q,T0,TolFun,maxiter,varargin)
% 
%  Args:
% 
%     - **Alead** [n x n x h x h array] : The suffix 01 is there to
%       indicate that the Aplus matrices are multiplied with the transition
%       probabilities: Jacobian of lead variables
% 
%      - **Acurr** [n x n x h array] : Jacobian of contemporaneous variables in
%        each regime
% 
%      - **Alag** [n x n x h array] : Jacobian of lagged variables in each
%        regime
% 
%      - **Q** [h x h matrix] : Transition matrix (used only if refinement)
% 
%      - **T0** [n x n x h array] : Initial guess for the solution (Not used)
% 
%      - **TolFun** [numeric] : Tolerance criterion for solution used in the
%        refinement of the solution
% 
%      - **maxiter** [numeric] : Maximum number of iterations used in the
%        refinement of the solution
% 
%      - **varargin** [] : additional arguments
%         - **refine** [true|{false}] : solve for sigma=1 instead of just
%           sigma=0
%         - **checkStab** [true|{false}] : check the stability of the
%            solution by the eigenvalues
%         - **allSols** [true|{false}] : flag for finding all solutions
%         - **msvOnly** [{true}|false] : return only the minimum state
%           variable solutions
%         - **xplosRoots** [{true}|false] : if all solutions are
%           computed we still can restrict ourselves to solutions that do not
%           involve explosive roots
%         - **debug** [true|{false}] : debug or not
% 
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
%  See also :  dsge_groebner, dsge_udc
%