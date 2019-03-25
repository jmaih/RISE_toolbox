%  diffuse_lyapunov_equation attempts to solve V=T*V*T'+RQR under
%  nonstationarity
% 
%  ::
% 
%    [P0,retcode]=diffuse_lyapunov_equation(T,RQR)
%    [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor)
%    [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor,check_unit_roots)
%    [P0,retcode]=diffuse_lyapunov_equation(T,RQR,diffuse_factor,check_unit_roots,diffuse_all)
% 
%  Args:
%     - T :
%     - RQR :
%     - diffuse_factor : [ positive scalar | {10} ]
%     - check_unit_roots : [ {true} | false ]
%     - diffuse_all : [ true | {false} ]
% 
%  Returns:
%     :
%     - P0 :
%     - Retcode :
% 
%  Note:
% 
%  Example:
% 
%     See also:
%