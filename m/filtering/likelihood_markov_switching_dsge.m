function [LogLik,Incr,retcode,obj]=likelihood_markov_switching_dsge(params,obj)%
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


if nargin~=2
    
    error([mfilename,':: Number of arguments must be 2'])
    
end

obj=assign_estimates(obj,params);

[obj,retcode]=solve(obj);

if ~retcode
    
    [obj,LogLik,Incr,retcode]=filter(obj);
    
    % evaluate the posterior
    if obj.options.debug
        
        disp(LogLik)
        
    end
    
end

if retcode
    
    Incr=[];
    
    LogLik=-obj.options.estim_penalty;
    
end

