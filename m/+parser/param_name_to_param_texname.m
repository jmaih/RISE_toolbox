function pname=param_name_to_param_texname(pname,chain_names)
% INTERNAL FUNCTION: change the parameter names from name_chain_state to name(chain,state)
%
% ::
%
%   pname=param_name_to_param_texname(pname,chain_names)
%
% Args:
%
%    - **pname** [char|cellstring: names of the parameters to change
%
%    - **chain_names** [cellstring]: names of the markov chains
%
% Returns:
%    :
%
%    - **pname** [char|cellstring]: names of the changed parameters
%
% See also:
%    - parser.param_texname_to_param_name
%

chain_names=cell2mat(strcat(chain_names(:)','|'));

pattern=['\_(',chain_names(1:end-1),')\_(\d+)'];

pname=regexprep(pname,pattern,'($1,$2)');

end