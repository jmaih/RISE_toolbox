function obj=setup_restrictions(obj)
% setup_restrictions - sets up linear, nonlinear and general restrictions
%
% Syntax
% -------
% ::
%
%	obj=setup_general_restrictions(obj)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: RISE object
%
% Outputs
% --------
%
% - **obj** [rise|dsge|rfvar|svar]: RISE object
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: RISE_GENERIC.SETUP_LINEAR_RESTRICTIONS,
% RISE_GENERIC.SETUP_NONLINEAR_RESTRICTIONS,
% RISE_GENERIC.SETUP_GENERAL_RESTRICTIONS 

if ~obj.restrictions_are_absorbed
    
    obj=setup_linear_restrictions(obj);
    
    obj=setup_nonlinear_restrictions(obj);
    
    obj=setup_general_restrictions(obj);
    
    obj.restrictions_are_absorbed=true;
    
end

end
