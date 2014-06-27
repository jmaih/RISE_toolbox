function [LogLik,Incr,retcode,obj]=likelihood_optimal_simple_rule(params,obj)

if nargin~=2
    error([mfilename,':: Number of arguments must be 2'])
end

if obj.NumberOfRegimes>1
    error([mfilename,':: OSR for markov switching not implemented yet'])
end
if obj.options.order>1
    error([mfilename,':: OSR with anticipations not implemented yet'])
end
%% this important output is not created yet
Incr=[];

%% solve the dsge model
[obj,retcode]=solve(obj,'evaluate_params',params);

%% initialize those and return if there is a problem
if retcode
    LogLik=-obj.options.Penalty;
else
    %% theoretical autocovariances. At order 0, we have variances only and we have already solved the model
    [A,retcode]=theoretical_autocovariances(obj,'ar',0);
    if retcode
        LogLik=-obj.options.Penalty;
    else
        % I could accelerate the step below by picking only the relevant
        % variables...
        weight_mat=diag(obj.W(:,:,1));
        if all(all(weight_mat==0))
            disp('Policy does not penalize deviations, you may need to check your loss function')
        end
        % loss is a positive number and we want to maximize the negative of
        % the loss. I do not have control over the way the user specifies
        % the loss function and so I first take the absolute value before
        % taking the negative(for maximization)
        Incr=-abs(weight_mat).*diag(A);
        % make sure we are going to do a minimization. We may not have
        % control over the way the user specified the objective function
        LogLik=sum(Incr);
        if isnan(LogLik)
            LogLik=-obj.options.Penalty;
        end % if isnan(loglik)
        if obj.options.kf_filtering_level
            error([mfilename,':: restrictions or filtering not implemented for optimal simple rules estimation'])
        end
    end % if retcode
end % if retcode



