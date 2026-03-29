%--- help for dsge/resimulate ---
%
%  Uses the intial conditions given in x0h and shocks information in
%  structure fx to resimulate the data using model m. The variables in m not
%  appearing in fx are initialized at zero if x0h is empty
% 
%  Syntax::
% 
%    mysimul=resimulate(m,fx)
% 
%    mysimul=resimulate(m,fx,x0h)
% 
%    [mysimul,retcode]=resimulate(...)
% 
%  Inputs:
% 
%  - m : [rise|dsge] model(s) to simulate
% 
%  - fx : (possibly modified output of) filtration.
% 
%  - x0h : (n_x x 1 x h): initial conditions in all the regimes. N.B: this
%    has to be in log-form. If the model has variables declared as
%    log_variables, those variables should be entered here in log. RISE will
%    not take the log of x0h.
% 
%  - fxlogvars : vector of booleans specifying which endogenous variables in
%    fx are log variables
% 
%  Outputs:
% 
%  - mysimul : [struct|cell array] structure or cell array of structures
%    with the simulations
% 
%  - retcode : [double] return code. 0 if there is no issue
% 
%  Note: all the simulations are returned in their log-form i.e. they are
%  not re-exponentiated !!!
% 
%  See also: dsge/historical_decomposition
%