function obj=log_marginal_data_density_mhm(obj)
if isempty(obj),obj=struct();return,end

simulation_folder=[obj.options.results_folder,filesep,'simulations'];
[theta_mean,V0]=parameters_posterior_moments(simulation_folder);
W = what(simulation_folder);
W=W.mat;
locs=find(strncmp('chain_',W,6));
if isempty(locs)
    error([mfilename,':: no simulations found'])
end
W=W(locs);
number_of_matrices=numel(W);
npar=size(theta_mean,1);
Vi=V0\eye(npar);
detV=det(V0);
tau=transpose(0.1:.1:.9);
MHM=0;
iter=0;
warn_state=warning;
warning off %#ok<WNOFF>
for m=1:number_of_matrices
    tmp=load([simulation_folder,filesep,W{m}]);
    Params=tmp.Params;
    log_post=-tmp.minus_logpost_params;
    nvals=size(Params,2);
    for ii=1:nvals
        lnf=geweke_truncated_multivariate_normal(Params(:,ii),theta_mean,detV,Vi,tau);
        MHM=MHM+exp(lnf-log_post(ii));
        iter=iter+1;
    end
end
warning(warn_state)
MHM=1./(MHM/iter);
log_mdd=log(MHM);
obj=obj.set_properties('log_mdd_mhm',max(log_mdd));

function lnf=geweke_truncated_multivariate_normal(theta,theta_mean,detV,Vi,tau)
d=numel(theta);
[lik,v_iF_v]=conditional_likelihood(theta-theta_mean,Vi,detV,d);
f=lik./tau.*(v_iF_v<=chi2inv(tau,d));
lnf=log(f);
