function [this,postSims]=posterior_simulator(obj,varargin)
% posterior_simulator -- simulate the posterior distribution of parameters
%
% Syntax
% -------
% ::
%
%   this=posterior_simulator(obj,varargin)
%
%   [this,postSims]=posterior_simulator(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: initial model object
%
% - **varargin** [pairwise addional inputs]: The most relevant for this
% function are:
%   - **mcmc_max_number_of_vectors_per_file** [integer|{1000}]:
%   - **mcmc_target_range** [interval|{[.25,.45]}]:
%   - **mcmc_burn_rate** [percentage|{10/100}]:
%   - **mcmc_number_of_simulations** [integer|{20000}]:
%   - **mcmc_initial_covariance_tune** [positive scalar|{0.25}]:
%   - **mcmc_number_of_parallel_chains** [integer|{2}]:
%   - **mcmc_delay_rejection** [false|{true}]:
%   - **mcmc_adapt_covariance** [true|{false}]:
%   - **mcmc_gibbs_sampler_func** [handle|{''}]:
%   - **mcmc_thin** [integer|{1}]: store every mcmc_thin draw
%   - **mcmc_diagcov_adjust_coef** [positive scalar|{1e-5}]:
%   - **mcmc_retune_cov_every** [integer|{100}]:
%
% Outputs
% --------
%
% - **this** [rise|dsge|svar|rfvar]: initial model object in which are
% added some posterior statistics. But not the computation of the marginal
% likelihood, which has to be done separately.
%
% - **postSims** [cell array]: place holder for posterior simulations.
%
% More About
% ------------
%
% - If the function is called only with one output, the draws are saved to
% disc. Else, they are returned in the second output argument.
%
% - uses a Metropolis Hastings algorithm with adaptation and delayed
% rejection. It is possible to start simulating the posterior distribution
% even without having maximized the posterior. 
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    this=struct('mcmc_max_number_of_vectors_per_file',1000,...
        'mcmc_target_range',[.25,.45],...
        'mcmc_burn_rate',10/100,...
        'mcmc_number_of_simulations',20000,... 
        'mcmc_initial_covariance_tune',0.25,...
        'mcmc_number_of_parallel_chains',2,...
        'mcmc_delay_rejection',true,...
        'mcmc_adapt_covariance',false,...
        'mcmc_gibbs_sampler_func','',...
        'mcmc_thin',1,... 10
        'mcmc_diagcov_adjust_coef',1e-5,... % COVt=COV_{t-1}+mcmc_diagcov_adjust_coef*eye(npar)
        'mcmc_retune_cov_every',100);
    % 'mcmc_diagnostics',false,...
    return
end

nn=length(varargin);
if rem(nn,2)
    error([mfilename,':: arguments must come in pairs'])
end

obj=set(obj,varargin{:});
% now make a copy so as to keep things inside intact
%----------------------------------------------------
this=obj;

obj=set(obj,'kf_filtering_level',0);
% the restrictions below may change the filtering level
%-------------------------------------------------------
obj=setup_linear_restrictions(obj);
obj=setup_general_restrictions(obj);
obj.estimation_under_way=true;

do_not_save_to_disk=nargout>1;

options=obj.options;
mcmc_number_of_simulations=options.mcmc_number_of_simulations;
mcmc_max_number_of_vectors_per_file=min(mcmc_number_of_simulations,...
    options.mcmc_max_number_of_vectors_per_file);
mcmc_number_of_parallel_chains=options.mcmc_number_of_parallel_chains;
mcmc_retune_cov_every=options.mcmc_retune_cov_every;

npar=numel(obj(1).estimation.priors);

number_of_matrices=ceil(mcmc_number_of_simulations/mcmc_max_number_of_vectors_per_file);
sampling_modes=[repmat({nan(npar,1)},1,mcmc_number_of_parallel_chains)
    repmat({nan},1,mcmc_number_of_parallel_chains)];
if do_not_save_to_disk
    postSims=cell(number_of_matrices,mcmc_number_of_parallel_chains);
end
post_sim_start_time=clock;

funevals=cell(1,mcmc_number_of_parallel_chains);
for pc=1:mcmc_number_of_parallel_chains
    [start,sampler,total_draws]=initialize_posterior_simulation(obj);
    utils.plot.waitbar('init',...
        struct('name',sprintf('MH for chain # %0.0f',pc),...
        'message','starting now...'));
    vcov_pc=0;
    xbar=0;
    matrix_index=0;
    sampling_modes{2,pc}=start.f0;
    sampling_modes{1,pc}=start.x0;
    iter_draw=0;
    Params=nan(npar,mcmc_max_number_of_vectors_per_file);
    minus_logpost_params=nan(1,mcmc_max_number_of_vectors_per_file);
    remain=mcmc_number_of_simulations;
    n_iter=[];
    offset=0;
    remain_in_one_run=mcmc_max_number_of_vectors_per_file;
    while remain
        nsamples=min([remain
            mcmc_max_number_of_vectors_per_file
            mcmc_retune_cov_every
            remain_in_one_run]);
        %------------------------------------------------------------------
        [Params(:,offset+(1:nsamples)), minus_logpost_params(offset+(1:nsamples)),acceptance_rate,start] = ...
            sampler(start,nsamples,...
            'waitbar_update',@(varargin)utils.plot.waitbar('update',percentage_update(),varargin{:}));
        %------------------------------------------------------------------
        
        % update the best
        %-----------------
        update_best(pc);
        
        % update the moments
        %-------------------
        [xbar,vcov_pc,n_iter]=utils.moments.recursive(xbar,vcov_pc,Params(:,offset+(1:nsamples)),n_iter);
        
        % update the cCs
        %---------------
        start=update_posterior_simulation_initial_conditions(obj,start,vcov_pc,acceptance_rate) ;
         
        % save down
        %----------
        offset=offset+nsamples;
        
        remain=max(0,remain-nsamples);
        remain_in_one_run=mcmc_max_number_of_vectors_per_file-offset;
        if remain_in_one_run==0 || remain == 0
            matrix_index=matrix_index+1;
            chain_number_matrix=sprintf('chain_%0.0f_%0.0f',pc,matrix_index);
            if do_not_save_to_disk
                postSims{matrix_index,pc}={chain_number_matrix;minus_logpost_params(1:offset);Params(:,1:offset)};
            else
                utils.parallel.par_save([simulation_folder,filesep,chain_number_matrix],...
                    {Params(:,1:offset),minus_logpost_params(1:offset)},{'Params','minus_logpost_params'})
            end
            offset=0;
            remain_in_one_run=mcmc_max_number_of_vectors_per_file;
        end
        
    end
    utils.plot.waitbar('close');
    funevals{pc}=start.funevals;
end

% obj.estimation_under_way=false; no longer required since the name of the
% output is changed

if do_not_save_to_disk
    tmp=postSims;
    postSims=struct();
    while size(tmp,2)
        col=tmp(:,1);
        tmp=tmp(:,2:end);
        for irow=1:numel(col)
            row_name=col{irow}{1};
            postSims.(row_name).minus_logpost_params=col{irow}{2};
            postSims.(row_name).Params=col{irow}{3};
        end
    end
    [theta_mean,theta_median,V0]=utils.mcmc.parameters_moments(postSims);
else
    [theta_mean,theta_median,V0]=utils.mcmc.parameters_moments(simulation_folder);
end

max_sim_id=cell2mat(sampling_modes(2,:));
max_sim_id=find(max_sim_id==max(max_sim_id),1,'first');
post_sim_mode=sampling_modes{1,max_sim_id};
post_sim_fmode=sampling_modes{2,max_sim_id};

post_sim_end_time=clock;

this.estimation.posterior_simulation=orderfields(...
    struct('mode',post_sim_mode,...
    'log_post',post_sim_fmode,'mean',theta_mean,'median',theta_median,'vcov',V0,...
    'mode_stdev',sqrt(diag(V0)),'log_marginal_data_density',struct(),...
    'post_sim_start_time',post_sim_start_time,'post_sim_end_time',post_sim_end_time,...
    'funevals',sum(cell2mat(funevals)))...
    );

    function update_best(pc)
        [~,tags]=sort(minus_logpost_params);
        if minus_logpost_params(tags(1))<sampling_modes{2,pc}
            sampling_modes{2,pc}=minus_logpost_params(tags(1));
            sampling_modes{1,pc}=Params(:,tags(1));
        end
    end

    function x=percentage_update(~)
        iter_draw=iter_draw+1;
        x=iter_draw/total_draws;
    end

end
