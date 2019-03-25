%  INTERNAL FUNCTION: Creates a handle to an m-file not on the matlab search path
% 
%  ::
% 
%    hdl = func2fhandle(full_path/func)
%    hdl = func2fhandle(relative_path/func)
%    [hdl1,hdl2,...,hdln] = func2fhandle(p1,p2,...,pn)
%    [hdl1,hdl2,...,hdln] = func2fhandle({p1,p2,...,pn})
% 
%  Args:
% 
%     - **pth** [char|cellstr] : path names to create handles for. The path
%       could be relative to the current directory or absolute
% 
%  Returns:
%     :
% 
%     - **hdl** [function_handle] : function handle(s)
% 
%