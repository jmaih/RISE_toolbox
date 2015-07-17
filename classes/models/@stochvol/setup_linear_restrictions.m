function obj=setup_linear_restrictions(obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


if isempty(obj)
    obj=setup_linear_restrictions@rise_generic(obj);
else
    error('linear restrictions on stochastic volatility model need to be updated')
    % for the stochastic volatility obj, I may still want to do as before, i.e.
    % applying the zero restrictions to determine the list of the estimated
    % parameters and not the other way around as it is done here, i.e. use the
    % list of estimated parameters to determine the restriction matrices.
end
    
end