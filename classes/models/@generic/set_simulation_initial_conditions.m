%--- help for set_simulation_initial_conditions ---
%
%  Set the initial conditions for forecasting, simulation and irfs
% 
%  ::
% 
%    Initcond=set_simulation_initial_conditions(obj)
% 
%  Args:
% 
%     obj (rise | dsge | svar | rfvar): model object
% 
%  Returns:
%     :
% 
%     - **Initcond** [struct]: Initial conditions for simulation
% 
%  Note:
% 
%     - if future values of endogenous are found, they are discarded the
%       variables are not explicitly declared as conditioning endogenous
%       variables.
% 
%     - if either endogenous or exogenous variables are declared as
%       conditioning variables, we automatically have a conditional forecasting
%       exercise.
% 
%     - if no variable is declared as conditioning, we have a simulation
%       exercise, in which case, there is a burn-in period.
% 
%     - Initialization is always done at the steady state and the declared
%       endogenous variables merely override the steady state values. This means
%       that if a variable is not found in the database, it is initialized at its
%       steady state.
% 
%  Example:
% 
%  See also:
%