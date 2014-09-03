function [obj,postSims]=posterior_simulator(obj,varargin)
% simulate the posterior distribution using a Metropolis Hastings algorithm
% with adaptation and delayed rejection. It is possible to start simulating
% the posterior distribution even without having maximized the posterior.
% the first output argument will be a RISE object containing some
% statistics from the posterior simulation. The computation of the marginal
% likelihood has to be done separately after this step.
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=struct('mcmc_max_number_of_vectors_per_file',1000,...
        'mcmc_target_range',[.25,.45],...
        'mcmc_burn_rate',10/100,...
        'mcmc_number_of_simulations',20000,... 55556 % <-- N*(1-burn_rate)/thin=5000
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
% obj=set(obj,'kf_filtering_level',0);

do_not_save_to_disk=nargout>1;

options=obj.options;
mcmc_number_of_simulations=options.mcmc_number_of_simulations;
mcmc_max_number_of_vectors_per_file=options.mcmc_max_number_of_vectors_per_file;
mcmc_burn_rate=options.mcmc_burn_rate;
mcmc_number_of_parallel_chains=options.mcmc_number_of_parallel_chains;
mcmc_thin=options.mcmc_thin;

number_of_burns=round(mcmc_burn_rate*mcmc_number_of_simulations);
% simulation_folder=obj.folders_paths.simulations;

npar=numel(obj(1).estimation.priors);

%----------------------------------------
lb=[obj.estimation.priors.lower_bound]';
ub=[obj.estimation.priors.upper_bound]';
drawfun=@(x,C)truncated_multivariate_normal.quick_and_dirty(x,C,lb,ub);
%----------------------------------------
number_of_matrices=ceil(mcmc_number_of_simulations/mcmc_max_number_of_vectors_per_file);
sampling_modes=[repmat({nan(npar,1)},1,mcmc_number_of_parallel_chains)
    repmat({nan},1,mcmc_number_of_parallel_chains)];
funevals=repmat({0},1,mcmc_number_of_parallel_chains);
if do_not_save_to_disk
    postSims=cell(number_of_matrices,mcmc_number_of_parallel_chains);
end
obj.estimation_under_way=true;
post_sim_start_time=clock;

[x0,f0,vcov0,objective,c0,adapt_covariance]=metropolis_initialize();
nsimu=number_of_burns+mcmc_number_of_simulations;
for pc=1:mcmc_number_of_parallel_chains
    utils.plot.waitbar('init',...
        struct('name',sprintf('MH for chain # %0.0f',pc),...
        'message','starting now...'));
    vcov_pc=vcov0;
    CS=transpose(chol(vcov_pc));
    c=c0;
    cCs=c*CS;
    fbest=f0;
    f_theta=f0;
    sampling_modes{2,pc}=f0;
    sampling_modes{1,pc}=x0;
    xbar=x0;
    theta=x0;
    Params=nan(npar,mcmc_max_number_of_vectors_per_file);
    minus_logpost_params=nan(1,mcmc_max_number_of_vectors_per_file);
    accepted=false(1,options.mcmc_retune_cov_every);
    retune_iter=0;
    store_index=0;
    matrix_index=0;
    ithin=0;
    newpeak=false;
    for d=1:nsimu
        [theta,f_theta]=metropolis_hastings_sampler(theta,f_theta);
        if d>number_of_burns
            ithin=ithin+1;
            if ithin==mcmc_thin||newpeak
                store_index=store_index+1;
                Params(:,store_index)=theta;
                minus_logpost_params(store_index)=f_theta;
                ithin=0;
                newpeak=false;
            end
            if options.mcmc_adapt_covariance
                [xbar,vcov_pc]=utils.moments.recursive(xbar,vcov_pc,theta,store_index+1);
            end
            if d==number_of_burns+mcmc_number_of_simulations
                % chop off and update t
                Params(:,store_index+1:end)=[];
                minus_logpost_params(store_index+1:end)=[];
                store_index=mcmc_max_number_of_vectors_per_file;
            end
            if store_index==mcmc_max_number_of_vectors_per_file
                % save down
                matrix_index=matrix_index+1;
                chain_number_matrix=sprintf('chain_%0.0f_%0.0f',pc,matrix_index);
                if do_not_save_to_disk
                    postSims{matrix_index,pc}={chain_number_matrix;minus_logpost_params;Params};
                else
                    utils.parallel.par_save([simulation_folder,filesep,chain_number_matrix],...
                        {Params,minus_logpost_params},{'Params','minus_logpost_params'})
                end
                % go to next
                store_index=0;
            end
        end
    end
    utils.plot.waitbar('close');
end

funevals=sum(cell2mat(funevals));
obj.estimation_under_way=false;

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

obj.estimation.posterior_simulation=orderfields(...
    struct('mode',post_sim_mode,...
    'log_post',post_sim_fmode,'mean',theta_mean,'median',theta_median,'vcov',V0,...
    'mode_stdev',sqrt(diag(V0)),'log_marginal_data_density',struct(),...
    'post_sim_start_time',post_sim_start_time,'post_sim_end_time',post_sim_end_time,...
    'funevals',funevals)...
    );

    function [theta,f_theta]=metropolis_hastings_sampler(theta0,f_theta0)
        retune_iter=retune_iter+1;
        [theta,f_theta,accepted(retune_iter),funevals{pc}]=...
            utils.mcmc.random_walk_mcmc(objective,theta0,f_theta0,drawfun,cCs,...
            options.mcmc_delay_rejection,funevals{pc});
        if -f_theta>-fbest
            newpeak=true;
            fprintf(1,'%s %5.0f %8.4f\n','new peak found in chain',pc,f_theta);
            fbest=f_theta;
            sampling_modes(:,pc)={theta;f_theta}';
        end
        
        if adapt_covariance
            [xbar,vcov_pc]=utils.moments.recursive(xbar,vcov_pc,theta,store_index+1);
        end
        if retune_iter==options.mcmc_retune_cov_every
            if adapt_covariance
                [R,p]=chol(vcov_pc+options.mcmc_diagcov_adjust_coef*eye(npar));
                if p==0
                    cCs=transpose(R);
                    fprintf(1,'%s %3.0d %s %8.4f %s \n',...
                        'chain',pc,'global peak',fbest,'covariance matrix adapted');
                end
            else
                [c,acceptance_rate]=utils.mcmc.retune_coefficient(c,options.mcmc_target_range,accepted);
                cCs=c*CS;
                fprintf(1,'%s %3.0d %s %8.4f %s %8.3f  %s %8.3f \n',...
                    'chain',pc,'global peak',fbest,'tunning coeff',c,'acceptance rate',acceptance_rate);
            end
            retune_iter=0;
            utils.plot.waitbar('update',d/nsimu);
        end
    end
    function [x0,f0,vcov,objective,c,adapt_covariance]=metropolis_initialize()
        c=obj.options.mcmc_initial_covariance_tune;
        if c<=0
            error([mfilename,':: mcmc_initial_covariance_tune must be positive'])
        end
        vcov=repair_covariance_matrix(obj.estimation.posterior_maximization.vcov);
        objective=@(x)fh_wrapper(x);
        x0=obj.estimation.posterior_maximization.mode;
        not_optimized=isempty(x0);
        if not_optimized
            x0=[obj.estimation.priors.start]';
            vcov=eye(npar);
            [minus_log_post,retcode]=objective(x0);
            f0=minus_log_post;
            if retcode
                msg=utils.error.decipher(retcode);
                error(['At MCMC simulation, starting values leads to the following problem ''',...
                    msg,...
                    ''' .This can be corrected. please contact junior.maih@gmail.com'])
            end
        else
            f0=-obj.estimation.posterior_maximization.log_post;
        end
        fprintf(1,'%s %8.4f\n','Initial value of log posterior ',f0);
        % adapt the covariance matrix automatically if we have a
        % diagonal original covariance matrix regardless of ...
        %--------------------------------------------------------------
        adapt_covariance=all(all(abs(diag(diag(vcov))-vcov)<1e-12));
        adapt_covariance=adapt_covariance||options.mcmc_adapt_covariance;
        function A_psd = repair_covariance_matrix(vcov)
            [V,D] = eig(vcov);
            A_psd = V*diag(max(diag(D),sqrt(eps)))*V';
        end
    end

    function [minus_log_post,retcode]=fh_wrapper(x)
        % this function returns the minimax if there are many objects
        % evaluated simultaneously
        
        retcode=0;
        if any(x<lb)||any(x>ub)
            retcode=7;
        else
            [~,minus_log_post,~,issue,viol]=...
                estimation_wrapper(obj,[],x,lb,ub,0); % function evaluations computed elsewhere
            if ~isempty(issue)
                retcode=issue;
            end
            if ~isempty(viol) && any(viol>0)
                retcode=7;
            end
        end
        if retcode
            minus_log_post=obj.options.estim_penalty;
        end
    end
end

