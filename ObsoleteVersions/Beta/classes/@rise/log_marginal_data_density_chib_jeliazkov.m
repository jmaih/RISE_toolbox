function obj=log_marginal_data_density_chib_jeliazkov(obj)
if isempty(obj),obj=struct();return,end

simulation_folder=obj.folders_paths.simulations;

% verbose=false;
x_mode=vertcat(obj.estimated_parameters.mode);
vcov_mode=obj.vcov;
npar=size(x_mode,1);
Vi=vcov_mode\eye(npar);
detV=det(vcov_mode);
W = what(simulation_folder);
W=W.mat;
locs=find(strncmp('chain_',W,6));
if isempty(locs)
    error([mfilename,':: no simulations found'])
end
W=W(locs);
number_of_matrices=numel(W);
CJ=0;
iter=0;
for m=1:number_of_matrices
    tmp=load([simulation_folder,filesep,W{m}]);
    Params=tmp.Params;
    nvals=size(Params,2);
    for ii=1:nvals
        lik=conditional_likelihood(x_mode-Params(:,ii),Vi,detV,npar);
        CJ=CJ+lik;
        iter=iter+1;
    end
end
CJ=CJ/iter;
% now sample J parameter vectors from the proposal density. I don't think I
% will save those draws... but perhaps I should put a seed to the random
% number generator such that I get the same results.
J=iter;
lb=vertcat(obj.estimated_parameters.lb);
ub=vertcat(obj.estimated_parameters.ub);

CS=transpose(chol(vcov_mode));
f0=obj.log_post;
alpha0=0;
for ii=1:J
    valid=false;
    while ~valid
        xi=x_mode+CS*randn(npar,1);
        if all(xi>=lb) && all(xi<=ub)
            valid=true;
        end
    end
    logpost=log_posterior_kernel(obj,xi);

	f_i=-logpost;
    a=alpha_probability(f_i,f0);
    alpha0=alpha0+a;
end
alpha0=alpha0/J;
log_mdd=log(CJ/alpha0);
obj=obj.set_properties('log_mdd_chib_jeliazkov',log_mdd);
