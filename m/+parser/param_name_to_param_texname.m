function pname=param_name_to_param_texname(pname,chain_names)
% PARAM_NAME_TO_PARAM_TEXNAME - change the parameter names from
% name_chain_state to name(chain,state)
%
% ::
%
%
%   pname=PARAM_NAME_TO_PARAM_TEXNAME(pname,chain_names)
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
% Note:
%
% Example:
%
%    See also: PARSER.PARAM_TEXNAME_TO_PARAM_NAME

chain_names=cell2mat(strcat(chain_names(:)','|'));

pattern=['\_(',chain_names(1:end-1),')\_(\d+)'];

pname=regexprep(pname,pattern,'($1,$2)');

end