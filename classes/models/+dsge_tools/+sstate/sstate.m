%--- help for dsge/sstate ---
%
%  Computes the steady state of a dsge model
% 
%  ::
% 
%    [obj,structural_matrices,retcode]=sstate(obj,varargin)
% 
%  Args:
% 
%     - obj (rise | dsge): model file
% 
%     - varargin : usual optional arguments in pairs
% 
%     - sstate_blocks (true|{false}): blockwise solution of the
%       steady state.
% 
%     - sstate_file (char | function_handle | {''}): name of the steady state file
% 
%     - sstate_use_file (false | {true}): use the steady state file or
%       program(in the model file) to solve the steady state 
% 
%     - sstate_solver (char | function_handle | {'lsqnonlin'}): Other solvers
%       known to RISE include "fsolve", "fminunc", "fminsearch". If you
%       want to supply your own solver, then it should be of the format
%       [x,f,exitflag]=yourSolver(Objective,x0,lb,ub,options,varargin) i.e.
%       it accepts the same inputs as lsqnonlin and returns the same 3
%       first outputs of fsolve
% 
%     - sstate_bgp_shift (numeric | {5}): shift/lead of the static model for
%       solving the balanced growth path.
% 
%     - sstate_unique (true | {false}): In a regime switching model,
%       there are potentially multiple steady states if some of the
%       switching parameters affect the steady state of the model. When
%       this is the case, the option "sstate_unique" is used for
%       forcing the steady state to be unique nevertheless. In that case,
%       the unique steady state is computed at the ergodic mean of the
%       parameters.
% 
%     - sstate_imposed (true | {false}): is typically used when 
% 
%       - one wants to approximate the model around a point that is not the
%         steady state of the model. For instance, in the piecewise linear
%         approach, there is a reference regime (normal times) which assumes
%         that no alternative regime exists and therefore the steady state of
%         the model is computed as the steady state of that particular
%         regime. That steady state is then imposed to be the steady state
%         even in the alternative regime (e.g. the zero-lower bound regime).
% 
%       - one wants to prevent RISE from treating the steady state file or
%         model as initial condition. For instance, under estimation, with
%         sstate_imposed set to false, RISE will always try to compute the
%         steady state numerically if the program computing the steady state
%         fails. But this is wasteful since the numerical solution can never
%         be better than the analytical solution. In other words, if the
%         analytical solution fails, the numerical solution cannot succeed.
% 
%       - one wants to prevent RISE from trying to solve the steady state if
%         the steady state equations do not hold. This useful for inspecting
%         the offending equations by means of the "resid" method. 
% 
%     - sstate_loop (true | {false}): This option is used to inform RISE
%       that the steady state file (or model) would accurately compute the
%       steady state of a selected number of variables if it is given the
%       steady state values for some other variables. Therefore, RISE will
%       loop over (iterate) over values of the unsolved variables until the
%       steady state is found.
%       This is typically the case for instance in optimal policy, where,
%       say, we do not know the steady state value of the interest rate.
%       Note that when sstate_loop is triggered, block decomposition
%       (sstate_blocks) becomes "inactive"
%       It is important to note that when the option is used in combination
%       with sstate_bounds, only the variables looped over are constrained.
% 
%     - sstate_use_jacobian (true | {false}): In linear models, the
%       jacobian is always used. In the nonlinear case, however, the true
%       jacobian tends not to work as well as its finite differences
%       approximation. Hence by default, we do not invoke the use of the true
%       jacobian.
% 
%     - sstate_endo_param_swap (cell array | {} ): When not empty, it is
%       a cell array with n rows and 3 columns, where n is the number of
%       restrictions, the first column gathers the names of the endogenous
%       variables whose values are given in the second column, the third column
%       includes the parameters names that are endogenized. 
% 
%     - sstate_impose_commitment (true|{false}): impose the steady state
%       of commitment to the solving of the steady state for discretion,
%       loose commitment and stochastic replanning.
% 
%     - sstate_bounds (struct | {} ): When not empty, it is
%       a structure such that each field is the name of an endogenous
%       variable for which we have constraints. Each field itself is either
% 
%       -  a scalar e.g bounds.C=2, which becomes bounds.C=[2,-inf,inf]
%       -  or a 3-element vector bounds.C=[2,-inf,inf], where the first
%          element represents the start value (initial guess), the second
%          element the lower bound and the third element the upper bound.
% 
%       It is also possible to simultaneously specify the growth rate of a
%       variable in case the model is nonstationary. In that case the
%       syntax above generalizes to a 1 x 2 cell array
%       bounds.C={sstateInfo,bgpInfo} where sstateInfo and bgpInfo at
%       scalars of 1 x 3 vectors as above.
% 
%       It is also possible to specify multiple regimes in case of a
%       regime-switching DSGE model. In that case the user as to specify
%       as many dimensions to the structure as the number of regimes.
%       e.g. bounds(2).C=bounds(1).C
% 
%       If no information is provided beyond the first regime, RISE
%       automatically extrapolates the first regime information to the
%       others.
% 
%       The default bounds for the steady states are -inf and inf
%       The default bounds for the growth rates are -log(2) and log(2),
%       i.e.::
% 
%          -log(2) <= C{t}-C{t-1} <= log(2) or
%          -log(2) <= log(C{t}/C{t-1}) <= log(2)
% 
%       Irrespective of the initial conditions chosen through sstate_bounds,
%       in the presence of a steady state file or a steady state model, the
%       procedure will go through such program and distort the initial
%       conditions
% 
%       Setting "sstate_bounds" to be initial_sstate(m), where m has been
%       previously solved, returns good initial conditions for both the
%       steady state and growth but poor bounds, which doesn't matter if
%       the guess is exact or good
% 
%       Setting "sstate_bounds" to be get(m,'sstate') returns the steady
%       state only, which is fine if the model is stationary, internally
%       those values are going to be imposed so that x0=lb=ub
% 
%       There is probably a need for a sstate_bgp or sstate_growth call that
%       will also return the growth rate simultaneously but this is not
%       implemented at this point
% 
%       In the presence of good or exact initial conditions, it is the
%       responsibility of the user to remove a defectuous or inexact steady
%       state program. This is done by simply reinitializing the steady
%       state m=set(m,'sstate_file','') This also suggests that having a
%       steady state model written inside the model file is not a good idea.
% 
%  Returns:
% 
%     - **obj** [rise|dsge]: model file
%     - **structural_matrices** [struct]: structure containing various
%       important elements for the solution of the system
%     - **retcode** [numeric]: 0 if there was problem computing the steady
%       state.
% 
%  Note:
% 
%     - There are 2 cases to consider:
% 
%       - The user does not provide any steady state equations: RISE will
%         attempt to solve the steady state using a vector of zeros as initial
%         guess. It makes life easy if the user provides the status of the
%         variables in the system i.e. whether they grow linear or log-linearly.
%       - The user provide some equations for solving the steady state. This is
%         done in two ways:
% 
%           1. the steady_state_model block: the variables that do not appear
%              in the block will be initialized at zero. Some parameters can
%              also be computed inside the block. The user can define an
%              optimization to solve for a subset of steady state values
%              inside the block. The block has three attributes:
% 
%               a. imposed(default=false): RISE computes the solution at the
%                  specified point without checking that the point solves for the
%                  steady state
%               b. unique (default=false): RISE computes the steady state at
%                  the ergodic distribution of the parameters. In case the
%                  probabilities are endogenous, the ergodic distribution of the
%                  parameters is itself a function of the steady state of the
%                  variables.
%               c. loop(default=false): RISE considers the equations
%                  calculating the steady state as true and just solves for the
%                  missing variables by looping over the steady state program. The
%                  user can then use the values pushed into the steady state
%                  program to calculate the steady state for the included
%                  variables.
% 
%           2. the steady state file: The user writes a function which can be
%              called in two possible ways:
% 
%                (i) [vnames,info]=ssfile();
% 
%                    In this case the first output argument is the list of variables
%                    for which the user computes the steady state; the second output
%                    is a structure with fields unique, imposed and loop
%                    just as in the case of the steady state model.
% 
%                (ii) The other call to the function is ::
% 
%                        [y,newp,retcode]=ssfile(y,p,d,id,obj)
% 
%                     In this case, the first input (y) is the vector of steady
%                     states, which is updated and returned as the first output. The
%                     locations of the modifications are indicated by the fourth
%                     input (id), which is computed based on the list of the
%                     variables in vnames above. As for the other outputs, p is a
%                     structure with parameters, d is a structure with definitions,
%                     obj is the model object in case the user needs some further
%                     information for computing the steady state. In case some
%                     parameters are computed in the steady state file, they should
%                     be returned in the structure "newp". The last output "retcode"
%                     indicates whether no problem was encoutered in the computation
%                     of the steady state (retcode=0) or the opposite (retcode =
%                     any number different from 0).
% 
%     - Writing the steady state file in this ways makes it possible to use the
%       same file whether there is regime switching or not.
% 
%  See also dsge.initial_sstate, dsge_tools.sstate.set_bounds
%