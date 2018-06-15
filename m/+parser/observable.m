function new=observable()
% INTERNAL FUNCTION
%

% format endogenous, parameters, observables, etc
new=struct('name',{},...
    'tex_name',{},...
    'is_endogenous',false,...
    'state_id',false);

