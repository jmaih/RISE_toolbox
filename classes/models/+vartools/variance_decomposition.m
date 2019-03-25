%--- help for abstvar/variance_decomposition ---
%
%  Compute the variance decompsition of the VAR
% 
%  ::
% 
%     [vd,retcode] = variance_decomposition(self,params,Rfunc,shock_names,varargin);
% 
%  Args:
%     self (var object):
%     params :(optional) parameter values
%     Rfunc (function handle): (optional) transform parameters into dynamics
%     shock_names : names of shocks
%     varargin
% 
%  Returns:
%     :
% 
%     - **vd** (struct): a struct containing the variance decomposition:
% 
%        - **infinity** (struct):
%        - **conditional** (struct):
% 
%     - **retcode** (integer): This passes through the retcode given by Rfunc. Currently (2018/07/19), this output always returns 0.
% 
%
%    Other functions named variance_decomposition
%
%       generic/variance_decomposition
%