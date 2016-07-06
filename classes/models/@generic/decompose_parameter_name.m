function [position,regime_states,pname,chain,state]=...
    decompose_parameter_name(obj,pname,~)
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

position=find(strcmp(pname,obj.parameters.name));

state=1;

chain='const';

% locate the state in the regimes
regime_states=1;

end
