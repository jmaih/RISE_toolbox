%--- help for dsge/indirect_inference ---
%
%  indirect_inference : Estimates coefficients using indirect inference
% 
%  ::
% 
%    [m]=indirect_inference(m,myobjective)
%    [m]=indirect_inference(m,myobjective,varargin)
% 
%  Args:
% 
%     m (rise | dsge): scalar or vector of model objects.
% 
%     myobjective (function handle): function that takes as input the model
%     object and returns :
% 
%        1. the criterion to minimize
%        2. a flag to inform whether the criterion is successfully computed
%            or not
% 
%     varargin : usual options for a dsge/rise object
% 
%  Returns:
%     :
% 
%     - **m** [scalar|vector]:  scalar or vector of model objects.
% 
%  Note:
%