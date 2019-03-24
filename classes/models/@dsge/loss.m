%--- help for dsge/loss ---
%
%  Calculates welfare
% 
%  ::
% 
%    [welf,retcode,V,d]=loss(obj)
%    [welf,retcode,V,d]=loss(obj,simuls)
%    [welf,retcode,V,d]=loss(obj,simuls,varargin)
% 
%  Args:
% 
%     obj (rise | dsge): scalar or vector of model objects.
% 
%     simuls (ts | struct | empty): if empty, the unconditional welfare is
%       returned. if not empty the conditional welfare is calculated instead
% 
%     varargin : optional arguments coming in pairs
% 
%  Returns:
%     :
% 
%     - **welf** [scalar|vector]: conditional loss (scalar) or unconditional
%       loss (scalar if number of regimes = 1)
% 
%     - **retcode** [scalar]: return code
% 
%     - **V** [n x n x h]: array of equilibrium matrices
% 
%     - **d** [h x 1 vector]: constant in loss
% 
%  Note:
% 
%     The loss is such that :math:`v(x_t,r_t)=x_t.'*Vrt*x_t+drt`
% 
%