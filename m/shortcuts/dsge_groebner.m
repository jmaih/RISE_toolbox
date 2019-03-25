%  dsge_groebner Groebner basis-based solution for DSGE models. Symbolic
%  solution of algebraic equations implied by a DSGE model. The procedure
%  attempts to find all possible solutions.
% 
%  ::
% 
%    [Tz_pb,eigvals,retcode]=dsge_groebner(Aplus01,A0,A_,Q,T0,TolFun,maxiter)  
% 
%  Args:
% 
%     - **Aplus01** [n x n x h x h array] : The suffix 01 is there to
%       indicate that the Aplus matrices are multiplied with the transition
%       probabilities: Jacobian of lead variables
% 
%      - **A0** [n x n x h array] : Jacobian of contemporaneous variables in
%        each regime
% 
%      - **A_** [n x n x h array] : Jacobian of lagged variables in each
%        regime 
% 
%      - **Q** [h x h matrix] : Transition matrix (Not used) 
% 
%      - **T0** [n x n x h array] : Initial guess for the solution (Not used) 
% 
%      - **TolFun** [numeric] : Tolerance criterion for solution (Not used) 
% 
%      - **maxiter** [numeric] : Maximum number of iterations (Not used) 
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
%  N:B:
%     :
% 
%     This solver is not guaranteed to give solutions and its computational
%     costs increase rapidly with the size of the problem. For instance I
%     have not been able to solve a problem with 18 unknowns. It is
%     therefore advisable to have a problem as small as possible.
% 
%     If Matlab shows busy for a long time, it is better just to kill it
%     rather than waiting.
% 
%     This procedures requires Matlab version 2018a or above.
% 
%  See also :  dsge_udc
%