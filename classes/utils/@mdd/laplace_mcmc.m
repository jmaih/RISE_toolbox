function log_mdd=laplace_mcmc(obj,~,~)

Hinv=cov(obj.theta_draws.');

log_mdd=utils.marginal_data_density.laplace_mdd(max(obj.LogPost_M),Hinv);

end
