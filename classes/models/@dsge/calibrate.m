%--- help for dsge/calibrate ---
%
%  calibrate : estimates the parameters of a model without the need of
%  evaluating the likelihood
% 
%  ::
% 
%    [m]=calibrate(m)
% 
%     m=calibrate(m,varargin)
% 
%  Args:
% 
%     - **m** (rise | dsge): model object or vectors of model objects
% 
%  Returns:
% 
%     - **m** (rise | dsge): model object or vectors of model objects with
%       information on the calibrated parameters
% 
%  .. index:: endogenous priors
% 
%  .. note::
% 
%     - It is assumed that varargin includes information on priors and
%       possibly endogenous priors
%