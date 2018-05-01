function [LogLik,Incr,retcode,obj]=likelihood_optimal_simple_rule(params,obj)
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

if obj.markov_chains.regimes_number>1
    
    error([mfilename,':: OSR for markov switching not implemented yet'])
    
end

if max(obj.exogenous.shock_horizon(:))>1
    
    error([mfilename,':: OSR with anticipations not implemented yet'])
    
end

% this important output is not created yet
%-----------------------------------------
Incr=[];

% solve the dsge model
%---------------------
obj=assign_estimates(obj,params);

[obj,retcode]=solve(obj);

if retcode
    
    LogLik=-obj.options.estim_penalty;
    
else
    
    [welf,retcode]=loss(obj);
    

    if retcode
        
        LogLik=-obj.options.estim_penalty;
        
    else
        
        % make sure we are going to do a minimization. We may not have
        % control over the way the user specified the objective function
        
        Incr=-abs(welf);
        
        LogLik=sum(Incr);
        
        if isnan(LogLik)
            
            LogLik=-obj.options.estim_penalty;
            
        end % if isnan(loglik)
        
        if obj.options.kf_filtering_level
            
            error([mfilename,':: restrictions or filtering not implemented for optimal simple rules estimation'])
        
        end
        
    end % if retcode

end % if retcode
