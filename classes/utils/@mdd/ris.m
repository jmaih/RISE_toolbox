function log_mdd=ris(obj,~,opts)

% Sylvia Frühwirth-Schnatter (2004): "Estimating marginal
% likelihoods for mixture and Markov switching models using bridge
% sampling techniques". Econometrics Journal 7(1) pp 143--167

Logq_M=old_draws_in_weighting_function(obj,opts,'RIS');

ratio = Logq_M-obj.LogPost_M;

ratiomax=max(ratio);

log_mdd = -1*(ratiomax+log(mean(exp(ratio-ratiomax))));

end
