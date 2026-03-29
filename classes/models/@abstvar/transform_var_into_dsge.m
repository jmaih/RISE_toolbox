%--- help for abstvar/transform_var_into_dsge ---
%
%  transform_var_into_dsge : transforms a parameterized VAR into a DSGE
%  model.
% 
%  SYNTAX
% 
%    m=transform_var_into_dsge(v)
% 
%    m=transform_var_into_dsge(v,exogenousValues)
% 
%  INPUTS
% 
%  - **v** [rfvar|svar|proxy_var|prfvar]: parameterized VAR object
% 
%  - **exogenousValues** [empty|struct]: structure containing the values for
%    the exogenous variables. This is because the DSGE does not contain
%    exogenous variables. Therefore the values of the exogenous variables
%    has to be fixed.
% 
%  OUTPUTS
% 
%  - **m** [rise]: parameterized DSGE model
% 
%  NB : At the moment, the routine does not support regime switching models.
%  The feature will be implemented if there is a demand for it.
%