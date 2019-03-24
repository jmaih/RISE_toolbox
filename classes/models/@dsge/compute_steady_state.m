%--- help for dsge/compute_steady_state ---
%
%  Computes the steady state of a dsge model
% 
%  ::
% 
%    [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
% 
%  Args:
% 
%     obj (rise | dsge): model file
% 
%     varargin : usual optional arguments
% 
%  Returns:
%     :
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
%               (a) imposed(default=false): RISE computes the solution at the
%                   specified point without checking that the point solves for the
%                   steady state
%               (b) unique (default=false): RISE computes the steady state at
%                   the ergodic distribution of the parameters. In case the
%                   probabilities are endogenous, the ergodic distribution of the
%                   parameters is itself a function of the steady state of the
%                   variables.
%               (c) loop(default=false): RISE considers the equations
%                   calculating the steady state as true and just solves for the
%                   missing variables by looping over the steady state program. The
%                   user can then use the values pushed into the steady state
%                   program to calculate the steady state for the included
%                   variables.
% 
%           2. the steady state file: The user writes a function which can be
%              called in two possible ways:
% 
%                (i) ::
% 
%                        [vnames,info]=ssfile();
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
%       same whether there is regime switching or not.
% 
%