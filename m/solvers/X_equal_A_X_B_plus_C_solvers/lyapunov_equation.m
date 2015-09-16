function [V,retcode]=lyapunov_equation(T,Q,options)
% lyapunov_equation solves the equation V=T*V*T'+Q
%
% Syntax
% -------
% ::
%   [V,retcode]=lyapunov_equation(T,Q)
%   [V,retcode]=lyapunov_equation(T,Q,options)
%
% Inputs
% -------
% - T :
% - Q :
% - options :
%
% Outputs
% --------
% - V :
% - retcode :
% 
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% example
% we would like to compute the variance of the system
% y_t=ay*y_{t-1}+ax*x_{t-1}+ae*e_t
% x_t=by*y_{t-1}+bx*x_{t-1}+be*e_t
% the matrix representation is 
% X_t=[ay ax;by bx]*X_{t-1}+[ae be]'*e_t;

% 
defaults=struct('lyapunov_algo','doubling',...
    'lyapunov_diffuse_factor',0);
defaults=utils.miscellaneous.mergestructures(doubling_solve(),defaults);
if nargin==0
	if nargout>1
		error([mfilename,':: number of output arguments cannot exceed 1 if there are no inputs'])
	end
	V=defaults;
	return
end

if nargin<3
    options=[];
end
if isfield(options,'lyapunov_algo')
    algo=options.lyapunov_algo;
else
    algo=defaults.lyapunov_algo;
end

switch algo
    case 'schur'
        error('Schur algorithms are under revision')
    case 'doubling'
        [V,retcode]=doubling_solve(T,[],Q,options);
    otherwise
        error([mfilename,':: unknown lyapunov algorithm option ',algo])
end

end

