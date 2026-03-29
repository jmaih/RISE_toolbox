%--- help for dsge/calculate_loss ---
%
%  Calculates welfare
% 
%  ::
% 
%    [welf,retcode,V,d]=calculate_loss(m,lossstr)
% 
%    [welf,retcode,V,d]=calculate_loss(m,lossstr,shocksdb)
% 
%    [welf,retcode,V,d]=calculate_loss(m,lossstr,shocksdb,varargin)
% 
%  Args:
% 
%     - **m** (rise | dsge): scalar or vector of model objects.
% 
%     - **lossstr** (char): loss function
% 
%     - **shocksdb** (ts | struct | empty): Shocks database to condition on over
%       simulation. if empty, the unconditional welfare is returned. if not
%       empty the conditional welfare is calculated instead
% 
%     - **varargin** : optional arguments coming in pairs
% 
%  Returns:
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