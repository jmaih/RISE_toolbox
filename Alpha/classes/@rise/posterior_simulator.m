function obj=posterior_simulator(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
obj=struct('mcmc_diagnostics',false,...
    'mcmc_max_number_of_vectors_per_file',1000,...
    'mcmc_target_range',[.25,.45],...
    'mcmc_burn_rate',10/100,...
    'mcmc_number_of_simulations',20000,...
    'mcmc_initial_covariance_tune',0.25,...
    'mcmc_number_of_parallel_chains',2);
    return
end

nn=length(varargin);
if rem(nn,2)
    error([mfilename,':: arguments must come in pairs'])
end

obj=set_options(obj,varargin{:});
obj=set_options(obj,'kf_filtering_level',0);

mcmc_initial_covariance_tune=obj.options.mcmc_initial_covariance_tune;
mcmc_target_range=obj.options.mcmc_target_range;
mcmc_number_of_simulations=obj.options.mcmc_number_of_simulations;
mcmc_max_number_of_vectors_per_file=obj.options.mcmc_max_number_of_vectors_per_file;
mcmc_diagnostics=obj.options.mcmc_diagnostics;
mcmc_burn_rate=obj.options.mcmc_burn_rate;
mcmc_number_of_parallel_chains=obj.options.mcmc_number_of_parallel_chains;

if mcmc_initial_covariance_tune<=0
	error([mfilename,':: mcmc_initial_covariance_tune must be positive'])
end
% hard-wired
returne_c_at=100;

number_of_burns=round(mcmc_burn_rate*mcmc_number_of_simulations);
simulation_folder=obj.folders_paths.simulations;

x0=obj.estimation.mode;
npar=size(x0,1);
lb=vertcat(obj.estimated_parameters.lb);
ub=vertcat(obj.estimated_parameters.ub);

CS=transpose(chol(obj.vcov));
number_of_matrices=ceil(mcmc_number_of_simulations/mcmc_max_number_of_vectors_per_file);
f0=obj.log_post;
sampling_modes=[repmat({x0},1,mcmc_number_of_parallel_chains)
    repmat({f0},1,mcmc_number_of_parallel_chains)];
parfor pc=1:mcmc_number_of_parallel_chains
    c=mcmc_initial_covariance_tune;
    cCS=c*CS;
    theta=x0;
    f_theta=f0;
    f00=f0;
    % burn initial draws
    accepted=false(1,returne_c_at);
    iter=0;
    for d=1:number_of_burns
        iter=iter+1;
        [theta,minus_logpost,accepted(iter)]=...
            random_walk_mcmc(obj,theta,f_theta,cCS,lb,ub);
        f_theta=-minus_logpost;
        if f_theta>f00
            fprintf(1,'%s %5.0f %8.4f\n','new peak found in chain',pc,f_theta);
            f00=f_theta;
            sampling_modes(:,pc)={theta;f_theta}';
        end
        if iter==returne_c_at
            [c,acceptance_rate]=retune_coefficient(c,mcmc_target_range,accepted);
            cCS=c*CS;
            iter=0;
            fprintf(1,'%s %3.0d %s %8.3f  %s %8.3f \n',...
                'chain(burn-in phase)',pc,'tunning coeff',c,'acceptance rate',acceptance_rate);
        end
    end
    % save the rest
    for m=1:number_of_matrices
        remains=mcmc_number_of_simulations-(m-1)*mcmc_max_number_of_vectors_per_file;
        n=min(mcmc_max_number_of_vectors_per_file,remains);
        Params=nan(npar,n);
        minus_logpost_params=nan(1,n);
        for d=1:n
            iter=iter+1;
            [Params(:,d),minus_logpost_params(d),accepted(iter)]=...
                random_walk_mcmc(obj,theta,f_theta,cCS,lb,ub);
            theta=Params(:,d);
            f_theta=-minus_logpost_params(d);
            if f_theta>f00
                fprintf(1,'%s %5.0f %8.4f\n','new peak found in chain',pc,f_theta);
                f00=f_theta;
                sampling_modes(:,pc)={theta;f_theta}'';
            end
            if iter==returne_c_at
                [c,acceptance_rate]=retune_coefficient(c,mcmc_target_range,accepted);
                cCS=c*CS;
                iter=0;
                fprintf(1,'%s %3.0d %s %8.4f %s %8.3f  %s %8.3f \n',...
                    'chain',pc,'global peak',f00,'tunning coeff',c,'acceptance rate',acceptance_rate);
            end
        end
        par_save([simulation_folder,filesep,'chain_',int2str(pc),'_',int2str(m)],...
            {Params,minus_logpost_params},{'Params','minus_logpost_params'})
    end
end

max_sim_id=cell2mat(sampling_modes(2,:));
max_sim_id=find(max_sim_id==max(max_sim_id),1,'first');
post_sim_mode=sampling_modes{1,max_sim_id};

% mean and quantiles
[theta_mean,~,quantiles]=parameters_posterior_moments(simulation_folder);
for ii=1:size(obj.estimated_parameters,1)    
    obj.estimated_parameters(ii)=obj.estimated_parameters(ii).set_properties(...
        'mean',theta_mean(ii),...
        'quantiles',quantiles{ii},...
        'post_sim_mode',post_sim_mode(ii));
end
if mcmc_diagnostics
    simulation_diagnostics(obj,simulation_folder);
end
end

function [x1,f1,accepted,alpha_prob]=random_walk_mcmc(obj,x0,f0,cCS,lb,ub)
npar=size(x0,1);
nobj=numel(obj);
theta_s=x0+cCS*randn(npar,1);
if any(theta_s<lb)||any(theta_s>ub)
    f_theta_s=-inf;
else
    minus_log_post=fh_wrapper(theta_s);
    f_theta_s=-minus_log_post;
end
alpha_prob=alpha_probability(f_theta_s,f0);
accepted=false;
if alpha_prob>rand
    x0=theta_s;
    f0=f_theta_s;
    accepted=true;
end
x1=x0;
f1=-f0;
        function [minus_log_post,retcode]=fh_wrapper(x)
        % this function returns the minimax if there are many objects
        % evaluated simultaneously
        
        fval=nan(1,nobj);
        % all objects are evaluated at the same point
        % the reason you want to output the object here is because it potentially
        % contains crucial information going forward. In particular, it contains
        % information about whether the model is stationary or not.
        for mo=1:nobj
            [fval(mo),~,~,~,retcode,obj(mo)]=log_posterior_kernel(obj(mo),x);
            if retcode
                break
            end
        end
        % Now take the negative for minimization
        minus_log_post=-min(fval);
    end
end

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