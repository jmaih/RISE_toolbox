function H = finite_differences(Objective,params,varargin)
% finite_differences - computes the hessian by finite differences
%
% ::
%
%
%	H = finite_differences(Objective,params)
%	H = finite_differences(Objective,params,varargin)
%
% Args:
%
%    - **Objective** [char|function handle]: function to differentiate
%
%    - **params** [vector]: point at which the differentiation is taken
%
%    - **varargin** : optional/further arguments of the objective function
%
% Returns:
%    :
%
%    - **H** [matrix]: Hessian matrix
%
% Note:
%
% Example:
%
%    See also:

H=utils.numdiff.hessian(Objective,params,[],varargin{:});