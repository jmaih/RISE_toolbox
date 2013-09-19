function [obj,postSims]=posterior_simulator(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=struct('mcmc_max_number_of_vectors_per_file',1000,...
        'mcmc_target_range',[.25,.45],...
        'mcmc_burn_rate',10/100,...
        'mcmc_number_of_simulations',20000,...
        'mcmc_initial_covariance_tune',0.25,...
        'mcmc_number_of_parallel_chains',2,...
        'mcmc_delay_rejection',true,...
        'mcmc_adapt_covariance',false,...
        'mcmc_diagcov_adjust_coef',1e-5,... % COVt=COV_{t-1}+mcmc_diagcov_adjust_coef*eye(npar)
        'mcmc_retune_cov_every',100);
% 'mcmc_diagnostics',false,...
            return
end

nn=length(varargin);
if rem(nn,2)
    error([mfilename,':: arguments must come in pairs'])
end

obj=set_options(obj,varargin{:});
obj=set_options(obj,'kf_filtering_level',0);

do_not_save_to_disk=nargout>1;

mcmc_initial_covariance_tune=obj.options.mcmc_initial_covariance_tune;
mcmc_target_range=obj.options.mcmc_target_range;
mcmc_number_of_simulations=obj.options.mcmc_number_of_simulations;
mcmc_max_number_of_vectors_per_file=obj.options.mcmc_max_number_of_vectors_per_file;
% mcmc_diagnostics=obj.options.mcmc_diagnostics;
mcmc_burn_rate=obj.options.mcmc_burn_rate;
mcmc_number_of_parallel_chains=obj.options.mcmc_number_of_parallel_chains;
mcmc_delay_rejection=obj.options.mcmc_delay_rejection;
mcmc_adapt_covariance=obj.options.mcmc_adapt_covariance;
mcmc_diagcov_adjust_coef=obj.options.mcmc_diagcov_adjust_coef;
mcmc_retune_cov_every=obj.options.mcmc_retune_cov_every;

if mcmc_initial_covariance_tune<=0
    error([mfilename,':: mcmc_initial_covariance_tune must be positive'])
end

number_of_burns=round(mcmc_burn_rate*mcmc_number_of_simulations);
simulation_folder=obj.folders_paths.simulations;

nobj=numel(obj);
x0=obj.estimation.posterior_maximization.mode;
vcov=obj.estimation.posterior_maximization.vcov;
npar=numel(obj.estimation.priors);
not_optimized=isempty(x0);
minus_log_post_func=@(x)fh_wrapper(x);
if not_optimized
    x0=[obj.estimation.priors.start]';
    vcov=eye(npar);
    % solve and establish stationarity
    %---------------------------------
    [minus_log_post,retcode]=minus_log_post_func(x0);
    f0=minus_log_post;
    if retcode
        msg=decipher_error(retcode);
        error(['At MCMC simulation, starting values leads to the following problem ''',...
            msg,...
            ''' .This can be corrected. please contact junior.maih@gmail.com'])
    end
else
    f0=-obj.estimation.posterior_maximization.log_post;
end
fprintf(1,'%s %8.4f\n','Initial value of log posterior ',f0);
% adapt the covariance matrix automatically if we have a diagonal original
% covariance matrix regardless of ...
%--------------------------------------------------------------------------
if all(all(abs(diag(diag(vcov))-vcov)<1e-12))
    mcmc_adapt_covariance=true;
end
lb=[obj.estimation.priors.lower_bound]';
ub=[obj.estimation.priors.upper_bound]';

number_of_matrices=ceil(mcmc_number_of_simulations/mcmc_max_number_of_vectors_per_file);
sampling_modes=[repmat({x0},1,mcmc_number_of_parallel_chains)
    repmat({f0},1,mcmc_number_of_parallel_chains)];
funevals=repmat({0},1,mcmc_number_of_parallel_chains);
if do_not_save_to_disk
postSims=cell(number_of_matrices,mcmc_number_of_parallel_chains);
end
post_sim_start_time=clock;
for pc=1:mcmc_number_of_parallel_chains
    vcov_pc=vcov;
    CS=transpose(chol(vcov_pc));
    c=mcmc_initial_covariance_tune;
    cCS=c*CS;
    theta=x0;
    f_theta=f0;
    f00=f0;
    xbar=x0;
    % burn initial draws
    %-------------------
    accepted=false(1,mcmc_retune_cov_every);
    retune_iter=0;
    for d=1:number_of_burns
        retune_iter=retune_iter+1;
        [theta,f_theta,accepted(retune_iter),funevals{pc}]=...
            random_walk_mcmc(minus_log_post_func,theta,f_theta,cCS,lb,ub,...
            mcmc_delay_rejection,funevals{pc});
        if -f_theta>-f00
            fprintf(1,'%s %5.0f %8.4f\n','new peak found in chain',pc,f_theta);
            f00=f_theta;
            sampling_modes(:,pc)={theta;f_theta}';
        end
        if retune_iter==mcmc_retune_cov_every
            [c,acceptance_rate]=retune_coefficient(c,mcmc_target_range,accepted);
            cCS=c*CS;
            retune_iter=0;
            fprintf(1,'%s %3.0d %s %8.3f  %s %8.3f \n',...
                'chain(burn-in phase)',pc,'tunning coeff',c,'acceptance rate',acceptance_rate);
        end
    end
    % save the rest
    %--------------
    for m=1:number_of_matrices
        remains=mcmc_number_of_simulations-(m-1)*mcmc_max_number_of_vectors_per_file;
        n=min(mcmc_max_number_of_vectors_per_file,remains);
        Params=nan(npar,n);
        minus_logpost_params=nan(1,n);
        for d=1:n
            retune_iter=retune_iter+1;
            [Params(:,d),minus_logpost_params(d),accepted(retune_iter),funevals{pc}]=...
                random_walk_mcmc(minus_log_post_func,theta,f_theta,cCS,lb,ub,...
                mcmc_delay_rejection,funevals{pc});
            theta=Params(:,d);
            f_theta=minus_logpost_params(d);
            if -f_theta>-f00
                fprintf(1,'%s %5.0f %8.4f\n','new peak found in chain',pc,f_theta);
                f00=f_theta;
                sampling_modes(:,pc)={theta;f_theta};
            end
            if mcmc_adapt_covariance
                t=(m-1)*mcmc_max_number_of_vectors_per_file+d;
                [xbar,vcov_pc]=rise_moments.recursive_moments(xbar,vcov_pc,Params(:,d),t+1);
            end
            if retune_iter==mcmc_retune_cov_every
                if mcmc_adapt_covariance
                    [R,p]=chol(vcov_pc+mcmc_diagcov_adjust_coef*eye(npar));
                    if p==0
                        cCS=transpose(R);
                        fprintf(1,'%s %3.0d %s %8.4f %s \n',...
                            'chain',pc,'global peak',f00,'covariance matrix adapted');
                    end
                else
                    [c,acceptance_rate]=retune_coefficient(c,mcmc_target_range,accepted);
                    cCS=c*CS;
                    fprintf(1,'%s %3.0d %s %8.4f %s %8.3f  %s %8.3f \n',...
                        'chain',pc,'global peak',f00,'tunning coeff',c,'acceptance rate',acceptance_rate);
                end
                retune_iter=0;
            end
        end
        chain_number_matrix=['chain_',int2str(pc),'_',int2str(m)];
        if do_not_save_to_disk
            postSims{m,pc}={chain_number_matrix;minus_logpost_params;Params};
        else
            par_save([simulation_folder,filesep,chain_number_matrix],...
                {Params,minus_logpost_params},{'Params','minus_logpost_params'})
        end
    end
end

funevals=sum(cell2mat(funevals));
if not_optimized
    funevals=funevals+1;
end
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
    [theta_mean,theta_median,V0]=parameters_posterior_moments(postSims);
else
    [theta_mean,theta_median,V0]=parameters_posterior_moments(simulation_folder);
end
% return
max_sim_id=cell2mat(sampling_modes(2,:));
max_sim_id=find(max_sim_id==max(max_sim_id),1,'first');
post_sim_mode=sampling_modes{1,max_sim_id};
post_sim_fmode=sampling_modes{2,max_sim_id};

% mcmc_diagnostics=true;
% if mcmc_diagnostics
%     simdiags=simulation_diagnostics(obj,simulation_folder);
% end
% 
post_sim_end_time=clock;

obj.estimation.posterior_simulation=orderfields(...
    struct('mode',post_sim_mode,...
    'log_post',post_sim_fmode,'mean',theta_mean,'median',theta_median,'vcov',V0,...
    'mode_stdev',sqrt(diag(V0)),'log_marginal_data_density',struct(),...
    'post_sim_start_time',post_sim_start_time,'post_sim_end_time',post_sim_end_time,...
    'funevals',funevals)...
    );

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

function [x1,f1,accepted,funevals,alpha_prob]=random_walk_mcmc(...
    minus_log_post_func,x0,f0,cCS,lb,ub,mcmc_delay_rejection,funevals)
npar=numel(x0);
[theta_s,minusLogPost_s]=new_proposal();
alpha_prob=alpha_probability(-minusLogPost_s,-f0);
accepted=alpha_prob>rand;
if accepted
    x1=theta_s;
    f1=minusLogPost_s;
else
    x1=x0;
    f1=f0;
    if mcmc_delay_rejection
        [theta_s2,minusLogPost_s2]=new_proposal();
        alpha_prob2=alpha_probability(-minusLogPost_s2,-minusLogPost_s);
        alpha13 = exp(-minusLogPost_s2-(-f0))*...
            qdens(theta_s2,theta_s)/qdens(x0,theta_s)*...
            (1-alpha_prob2)/(1-alpha_prob);
        accepted=alpha13>rand;
        if accepted
            x1=theta_s2;
            f1=minusLogPost_s2;
        end
    end
end
    function qab=qdens(a,b)
        ab=a-b;
        ab=ab(:);
        qab=exp(-0.5*(ab'*ab));
    end
    function [d,minusLogPost]=new_proposal()
        d=x0+cCS*randn(npar,1);
        if any(d<lb)||any(d>ub)
            minusLogPost=inf;
        else
            minusLogPost=minus_log_post_func(d);
            funevals=funevals+1;
        end
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