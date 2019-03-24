%  INTERNAL FUNCTION: Computes moments Ekron(u,u,...,u)
% 
%  ::
% 
%    [M1,M2,M3,...,M5]=u_higher_order_moments(siz)
% 
%  Args:
% 
%     - siz (struct): struct with fields
% 
%       - np: number of predetermined
%       - nb: number of predetermined and forward-looking
%       - ne: number of shocks
%       - nz: (total) number of state variables (pred,bobth,sig,shocks)
% 
%  Returns:
%     :
%     varargout
% 
%