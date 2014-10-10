function pname=valid_param_name_to_tex_name(pname,chain_names)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


% change the parameter names from name_chain_state to name(chain,state) and
% is therefore the opposite of param_name_to_valid_param_name

chain_names=cell2mat(strcat(chain_names(:)','|'));
pattern=['\_(',chain_names(1:end-1),')\_(\d+)'];
pname=regexprep(pname,pattern,'($1,$2)');