%  INTERNAL FUNCTION: change the parameter names from name(chain,state) to name_chain_state
% 
%  ::
% 
%    [pname_out,capture_errors]=PARAM_TEXNAME_TO_PARAM_NAME(pname)
% 
%  Args:
% 
%     - **pname** [char|cellstring]: names of the parameter names to change
% 
%  Returns:
%     :
% 
%     - **pname_out** [char|cellstring]: names of the changed parameter names
%     - **capture_errors** [cellstring]: list of invalid parameter names. If
%       this output is not requested, an error is issued for the very first
%       offending parameter name.
% 
%  See also:
%     parser.param_name_to_param_texname
% 
%