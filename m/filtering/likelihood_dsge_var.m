function [LogLik,Incr,retcode,obj]=likelihood_dsge_var(params,obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% N.B: This function assumes the likelihood is to be maximized

if nargin~=2
    error([mfilename,':: Number of arguments must be 0 or 2'])
end

% this important output is never created
Incr=[];

obj=assign_estimates(obj,params);

[obj,retcode,LogLik]=bvar_dsge(obj);

end 

