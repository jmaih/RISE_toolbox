function resid=hyperparameter_residuals(x,cdfn,plb,pub,prob,callingfun,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


a=x(1);
b=x(2);
alpha=1-prob;
resid=[cdfn(plb,a,b,varargin{:})-.5*alpha;
    cdfn(pub,a,b,varargin{:})-(1-.5*alpha)]; % matlab_gamma_definition
if any(isnan(resid))
    resid=1e+8*ones(2,1);
end
if ismember(callingfun,{'fmincon','fminsearch'})
    resid=norm(resid);
end
end
