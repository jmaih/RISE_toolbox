function vout=valid_names_in_text(vin)
% INTERNAL FUNCTION
%

% change the parameter names from name(chain,state) to name_chain_state
% from a text and not just from a list as in param_texname_to_param_name

vout=regexprep(vin,'(\w+)(\()([a-zA-Z]\w*),(\d+)(\))','$1_$3_$4');

