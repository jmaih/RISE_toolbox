%  INTERNAL FUNCTION: discards the options not used by a function
% 
%  ::
% 
%    opt=reselect_options(options,fn)
% 
%  Args:
% 
%     - **options** [struct]: all options
%     - **fn** [char|function_handle]: function whose options we are interested
%       in.
% 
%  Returns:
%     :
% 
%     - **opt** [struct]: options of function fn
% 
%  Note:
% 
%     - the function of interest should be such that if called without inputs,
%       it returns its default options.
% 
%