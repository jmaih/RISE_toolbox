function H = finite_differences(varargin)

H=utils.numdiff.hessian(varargin{:});