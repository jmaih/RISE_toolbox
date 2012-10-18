function [LogLik,Incr,retcode,obj]=likelihood_markov_switching_dsge(params,obj)%

if nargin~=2
    error([mfilename,':: Number of arguments must be 2'])
end

[obj,LogLik,Incr,retcode]=filter(obj,'evaluate_params',params);
% evaluate the posterior
if obj.options.debug
    disp(LogLik)
end
if retcode
    LogLik=-obj.options.Penalty;
end

