%  INTERNAL FUNCTION: flags the parameters that are modified in the steady state model
% 
%  ::
% 
%    is_changed=parameters_changed_in_ssmodel(shadow_sstate,pstring,n)
% 
%  Args:
% 
%     - **shadow_sstate** [cellstr]: list of steady state equations
%     - **pstring** [char]: generic name of shadow parameters
%     - **n** [integer]: number of parameters
% 
%  Returns:
%     :
% 
%     - **is_changed** [logical]: 1 x n vector flagging the changed
%       parameters
% 
%