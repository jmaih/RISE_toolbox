function vout=valid_names_in_text(vin)
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

% change the parameter names from name(chain,state) to name_chain_state
% from a text and not just from a list as in param_name_to_valid_param_name

vout=regexprep(vin,'(\w+)(\()([a-zA-Z]\w*),(\d+)(\))','$1_$3_$4');

