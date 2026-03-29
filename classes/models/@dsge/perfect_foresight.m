%--- help for dsge/perfect_foresight ---
%
%  perfect_foresight : Perfect foresight and extended path simulations
% 
%  ::
% 
% 
%    [db,states,retcode] = perfect_foresight(obj,varargin)
% 
%  Args:
% 
%     - **obj** [dsge|rise]: model object
% 
%     - **varargin** : additional arguments including but not restricted to
% 
%       - **simul_exogenous_func** [function handle|{empty}]: function that
%         generates exogenous values for exogenous variables
% 
%       - **simul_stack_solve_algo** [{'sparse'}|'fsolve'|'lsqnonlin'| user defined]:
%         alorithm for solving the problem. If the algorithm is user defined,
%         it should accept the same inputs as lsqnonlin and return the same
%         outputs as fsolve.
%         [y1,fval,exitflag,output]=user_algo(objfun,y0,lb,ub,options); 
%         output is a structure whose only relevant field is "iterations"
% 
%       - **simul_recursive** [true|{false}]: if true, extended path is
%         triggered. It is assumed that agents no longer anticipate the whole
%         sequence of shocks. Rather, in each period they expect zero shocks
%         in the future but are suprised that the realizations of the shocks
%         are nonzero. If false, perfect foresight is computed.
% 
%       - **simul_enforce_hybrid** [true|{false}]: if true, all subproblems
%         are solved assuming that there are both backward and forward-looking
%         terms. Use this option to avoid a long sequence of for-loops and
%         simple iterations in problems that can otherwise be solved one
%         period at a time. This happens for instance when the simulation
%         horizon is very long.
% 
%       - **simul_pf_bounds** [struct|{empty}]: sets the search space/bounds
%         for endogenous variables during perfect foresight and extended path
%         simulations. Those constraints can be implemented in similar fashion
%         to the way the constraints for the steady state are implemented but
%         using a distinct option rather than 'sstate_bounds'.
%         e.g. Assume we want to impose S>0, R<1 and I<1. This is implemented
%         as 
%         bounds=struct(); bounds.S=[nan,0,inf]; bounds.R=[nan,-inf,1]; bounds.I=[nan,-inf,1];
%         sims=perfect_foresight(...,'simul_pf_bounds',bounds)
% 
% 
%  Returns:
% 
%     - **sims** [struct]: simulated series
% 
%     - **fval** [scalar|vector]: maximum error
% 
%     - **retcode** [integer]: if 0, the simulation went fine. Else something
%       got wrong. In that case one can understand the problem by running
%       decipher(retcode)
% 
%  Note:
% 
%     - the algorithm allows for regime switching : the economy then jumps
%       from one regime to the other and the economic agents expect and are
%       aware of those jumps in advance. In this case the user should specify
%       the whole history of regimes
% 
%     - the algorithm allows for the use of symbolic, automatic or
%       numerical derivatives for the computation of the jacobian.
% 
%     The first element in each vector is supposed to represent the start
%     value but this is not used in this context. In a future update, those
%     start values will override whatever start values (or more precisely
%     the initial conditions) specified in the historical database
% 
%     - the algorithm uses block decomposition to efficiently solve the
%       simulation problem recursively. Besides the obvious advantage that
%       solving many small problems is more efficient than solving one large
%       one, it turns out that the decomposition will also isolate blocks that
%       are either purely backward-looking, purely forward-looking or purely
%       static. In those three cases further computational gains can be
%       obtained because then the problem can be solved period by period.
% 
%     - Perfect foresight is kind of nice in the sense that it illustrates
%       the fact that although the true solution of the model is not analytic,
%       it is potentially of the form x{t}=T(x{t-1},e{t},e{t+1},...,e{n}),
%       where n is the last observation. It is interesting to note that going
%       forward, we have x{t+1}=T(x{t},e{t+1},...,e{n}) ... so that in the end
%       x{n}=T(x{n-1},e{n})
% 
%     - Another insight from perfect foresight relates to the equivalence
%       between a stochastic model and a deterministic model when both of them
%       are linear : we get certainty equivalence, which can be verified in
%       the fact that the perfect foresight and the stochastic solutions yield
%       the same impulse response functions. Of course this happens only when
%       there is a unique solution. Under multiple stable solutions, there is
%       still certainty equivalence but the solution found using the perfect
%       foresight algorithm might depend on the vagaries of the optimization
%       routine.
% 
%  Example:
% 
%     See also: dsge.simulate
%