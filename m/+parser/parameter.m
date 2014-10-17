function new=parameter()
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
    'is_switching',false,...
    'is_in_use',false,...
    'governing_chain',{});
