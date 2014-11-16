function H = finite_differences(Objective,params,varargin)
% finite_differences - computes the hessian by finite differences
%
% Syntax
% -------
% ::
%
%	H = finite_differences(Objective,params)
%	H = finite_differences(Objective,params,varargin)
%
% Inputs
% -------
%
% - **Objective** [char|function handle]: function to differentiate
%
% - **params** [vector]: point at which the differentiation is taken
%
% - **varargin** : optional/further arguments of the objective function
%
% Outputs
% --------
%
% - **H** [matrix]: Hessian matrix
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

H=utils.numdiff.hessian(Objective,params,[],varargin{:});