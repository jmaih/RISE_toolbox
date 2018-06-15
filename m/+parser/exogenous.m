function new=exogenous()
% INTERNAL FUNCTION
%

% format endogenous, parameters, observables, etc
new=struct('name',{},...
    'tex_name',{},...
    'is_observed',false,...
    'is_in_use',false);

