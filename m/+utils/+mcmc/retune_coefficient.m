function [c,acceptance_rate]=retune_coefficient(mcmc_initial_covariance_tune,mcmc_target_range,accepted)
acceptance_rate=sum(accepted)/numel(accepted);
c=mcmc_initial_covariance_tune;
if acceptance_rate<mcmc_target_range(1)||acceptance_rate>mcmc_target_range(2)
    % increase c if the acceptance rate is high
    % decrease it otherwise
    target_rate=mean(mcmc_target_range);
    c=mcmc_initial_covariance_tune*acceptance_rate/target_rate;
end
end