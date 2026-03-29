%--- help for dsge/optimal_simple_rule ---
%
%  optimal_simple_rule : Estimates optimal simple rules coefficients
% 
%  ::
% 
%    [m]=optimal_simple_rule(m,lossstr)
%    [m]=optimal_simple_rule(m,lossstr,shocks)
%    [m]=optimal_simple_rule(m,lossstr,[],varargin)
%    [m]=optimal_simple_rule(m,lossstr,shocks,varargin)
% 
%  Args:
% 
%     m (rise | dsge): scalar or vector of model objects.
% 
%     lossstr (string|cell array): loss function two possibilities
%            - string : lossString. In this case discount = 0.99
%            - cell array : {discount,lossString}
% 
%     priors (struct): parameters to estimate
% 
%     shocks (ts | struct | double| function handle| empty): if empty, the
%       theoretical (unconditional) welfare will be calculated. Otherwise
%       the conditional welfare will be computed.
% 
%  Returns:
%     :
% 
%     - **m** [scalar|vector]:  scalar or vector of model objects.
% 
%  Note:
%