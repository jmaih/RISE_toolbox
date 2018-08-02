function obj=setup_restrictions(obj)
% INTERNAL FUNCTION: sets up linear, nonlinear and general restrictions
%
% ::
%
%	obj=setup_general_restrictions(obj)
%
% Args:
%
%    obj (rise | dsge | rfvar | svar): RISE object
%
% Returns:
%    :
%
%    - **obj** [rise|dsge|rfvar|svar]: RISE object
%
% See also:
%    - rise_generic.setup_linear_restrictions,
%    - rise_generic.setup_nonlinear_restrictions,
%    - rise_generic.setup_general_restrictions
%

if ~obj.restrictions_are_absorbed

    obj=setup_linear_restrictions(obj);

    obj=setup_nonlinear_restrictions(obj);

    obj=setup_general_restrictions(obj);

    obj.restrictions_are_absorbed=true;

end

end
