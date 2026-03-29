%--- order.m not found. Showing help for ode instead. ---
%
% ODE   Ordinary Differential Equation
%    F = ODE creates an ode object with default properties.
% 
%    F = ODE(Name=Value) sets the property NAME to VALUE.
% 
%    ODE properties:
%      Problem definition properties:
%                    ODEFcn - Function that defines the
%                             equations to solve: dy/dt = ODEFcn(t,y)
%               InitialTime - Initial time for integration. The default
%                             value is 0.
%              InitialValue - Value of the solution y at the initial time
%                Parameters - Parameters in the equations to pass
%                             to ODEFcn and other functions
%                  Jacobian - odeJacobian object defining the
%                             Jacobian matrix for ODEFcn
%           EventDefinition - odeEvent object defining events
%                             (zero-crossings) for the solver to detect
%                MassMatrix - odeMassMatrix object defining the mass
%                             matrix M to solve M*dy/dt = ODEFcn(t,y)
%      NonNegativeVariables - Indices of solution components
%                             that must be non-negative
%              InitialSlope - Value of dy/dt at the initial
%                             time, used by some solvers
% 
%      Solver control properties:
%                    Solver - Name of a solver, or "auto", "stiff", or
%                             "nonstiff" to automatically select a solver
%            SelectedSolver - Solver chosen when the solver is automatically
%                             selected with "auto", "stiff", or "nonstiff".
%                             This is read only. To manually update the
%                             solver, use the Solver property.
%             SolverOptions - Object that holds options
%                             specific to the selected solver
%         AbsoluteTolerance - Absolute error tolerance
%         RelativeTolerance - Relative error tolerance
% 
%      ODE methods:
%              solve - Solve the ODE over an interval or at specified points.
%        solutionFcn - Create function to interpolate the solution
%                      of the ODE object over an interval.
% 
%      Example: Sine and Cosine integration
%          D = ode;
%          % ODEFcn must return a column vector.
%          D.ODEFcn = @(t,y) [y(2); -y(1)];
%          D.InitialTime = 0;
%          % Initial values can be row or column vector
%          D.InitialValue = [0; 1];
%          sol = solve(D,-pi,pi)
%          sineFcn = solutionFcn(D,-pi,pi,OutputVariables=1)
%          tt = linspace(-pi,pi,19);
%          plot(sol.Time,sol.Solution,'-',tt,sineFcn(tt),'o');
% 
%       See also ode23, ode23s, ode23t, ode23tb, ode45, ode78, ode89,
%       ode113, ode15s
%
%    Documentation for ode
%       doc ode
%
%