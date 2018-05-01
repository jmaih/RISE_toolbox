function [log_post,log_lik,log_prior,Incr,retcode,obj,x1]=conclude_estimation(obj,x1)
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

if ~obj.is_optimal_simple_rule_model
    obj.options.kf_filtering_level=3;
end
[log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj,x1);

if obj.is_dsge_var_model
    obj.dsge_var.log_post=log_post(2);
    log_post=log_post(1);
    obj.dsge_var.log_lik=log_lik(2);
    log_lik=log_lik(1);
    obj.dsge_var.log_prior=log_prior;
end
if obj.options.kf_filtering_level
    % put the filters in the time series format
    obj=save_filters(obj);
end
% now we change the flag so that stability can be tested
% and parameters can be written back to their object
obj.estimation_under_way=false;

end