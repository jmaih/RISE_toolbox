function lmdd=log_marginal_data_density(obj,type,simulation_folder)
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

if isempty(obj)
lmdd=struct();
return
end

if nargin<3
    simulation_folder=[];
end
if isempty(simulation_folder)
    simulation_folder=obj.folders_paths.simulations;%[obj.options.results_folder,filesep,'simulations'];
end
is_saved_to_disk=ischar(simulation_folder);
if is_saved_to_disk
    W = what(simulation_folder);
    W=W.mat;
    locs=find(strncmp('chain_',W,6));
    if isempty(locs)
        error([mfilename,':: no simulations found'])
    end
    W=W(locs);
elseif isstruct(simulation_folder)
    W=fieldnames(simulation_folder);
else
    error('wrong specification of input')
end

% first determine the number of chains and the number  of matrices in each
% chain.
number_of_parallel_chains=obj.options.mcmc_number_of_parallel_chains;
switch type
    case 'mhm'
        xmean=obj.estimation.posterior_simulation.mean;
        vcov_mean=obj.estimation.posterior_simulation.vcov;
        self=modified_harmonic_mean(xmean,vcov_mean);
    case 'chib_jeliazkov'
        xmode=obj.estimation.posterior_simulation.mode;
        % below, here I am using the covariance matrix of the mean although
        % perhaps I should be using the covariance matrix of the mode
        vcov_mode=obj.estimation.posterior_simulation.vcov;
        log_post_mode=obj.estimation.posterior_simulation.log_post;
        lower_bound=obj.estimation.priors.lower_bound;
        upper_bound=obj.estimation.priors.upper_bound;
        log_post_func=@(x)log_posterior_kernel(obj,x);
        self=chib_jeliazkov(xmode,vcov_mode,log_post_mode,...
            lower_bound,upper_bound,log_post_func);
    otherwise
        error(['unknown marginal likelihood computation method ''',type,''])
end

for pc=1:number_of_parallel_chains
    matrices=regexp(W,['(?<!\w)chain_',sprintf('%0.0f',pc),'_\d+(?!\w)'],'match');
    matrices=[matrices{:}];
    number_of_matrices=numel(matrices);
    for m=1:number_of_matrices
        this_matrix=['chain_',sprintf('%0.0f',pc),'_',sprintf('%0.0f',m)];
        if is_saved_to_disk
            tmp=load([simulation_folder,filesep,this_matrix]);
        else
            tmp=simulation_folder.(this_matrix);
        end
        Params=tmp.Params;
        log_post=-tmp.minus_logpost_params;
        nvals=size(Params,2);
        for ii=1:nvals
            self=update(self,Params(:,ii),log_post(ii));
        end
    end
end

lmdd=conclude(self);


