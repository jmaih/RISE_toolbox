function [log_post,log_lik,log_prior,Incr,retcode,obj]=conclude_estimation(obj,x1)

obj.estimation_under_way=false;

[log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj,x1);

end