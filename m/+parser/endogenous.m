function new=endogenous()
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


% format endogenous, parameters, observables, etc
new=struct('name',{},...
    'tex_name',{},...
    'is_original',true,...
    'is_lagrange_multiplier',false,...
    'is_static',false,...
    'is_predetermined',false,...
    'is_pred_frwrd_looking',false,...
    'is_frwrd_looking',false,...
    'is_log_var',false,...
	'is_auxiliary',false);
