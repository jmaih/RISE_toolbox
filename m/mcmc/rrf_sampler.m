function [Results]=rrf_sampler(logf,lb,ub,options,mu,SIG)
% RRF_SAMPLER -- Rapid Reaction Force Sampler
%
% ::
%
%
%   [Results]=RRF_SAMPLER(logf,lb,ub)
%
%   [Results]=RRF_SAMPLER(logf,lb,ub,options)
%
%   [Results]=RRF_SAMPLER(logf,lb,ub,options,mu)
%
%   [Results]=RRF_SAMPLER(logf,lb,ub,options,mu,SIG)
%
% Args:
%
%    - **logf** [char|function_handle]: Objective function to MINIMIZE!!!
%
%    - **lb** [d x 1 vector]: lower bound of the paramters
%
%    - **ub** [d x 1 vector]: upper bound of the paramters
%
%    - **options** [struct]:
%      - **c** [scalar|{1}]: initial scale for the covariance matrix
%      - **c_range** [vector|{[sqrt(eps),100]}]: range of variation of c
%      - **alpha** [scalar|2-element|{[.25,.45]}]: target acceptance rate
%      - **burnin** [integer|{0}]: number of burn-in initial simulations
%      - **N** [integer|{20000}]: number of simulations
%      - **M_over_N** [numeric|{5/100}]: percentage of draws per bin
%      - **H** [integer|{40}]: number of tempering stages
%      - **extra_runs** [integer|{5}]: number of extra runs after lambda=1
%      - **lambda_1** [numeric|{2.5e-5}]: initial tempering
%      - **ps** [numeric|{0.9}]: probability of storing a draw
%      - **p** [numeric|{[]}]: probability of drawing from previous
%      distribution
%      - **p_mutant** [numeric|{0}]: probability of drawing a mutant
%      - **ess_min** [numeric|{0.1}]: constant for controlling the evolution
%      of lambda
%      - **focus_fire_after** [numeric|{0.5}]: percentage after which
%      exploitation takes over exploration
%      - **penalty** [numeric|{1e+8}]: worst possible function value
%      - **fixed_scaling** [true|{false}]: if true, the scaling (c) of the
%      covariance matrix is kept constant
%      - **geometric_lambda** [{true}|false]: geometric evolution of lambda
%      - **rwm_exp** [numeric|{0.6}]: tuning hyper-parameter for scale and
%      covariance matrix
%      - **use_true_moments** [true|{false}]: if true, the updated exact
%      covariance matrix of the draws is used at each step. If false, a
%      different update of the covariance matrix is used.
%
%    - **mu** [d x 1 vector]: initial condition for the sampler
%
%    - **SIG** [d x d matrix]: initial covariance matrix
%
% Returns:
%    :
%
%    - **Results** [struct]:
%      - **pop** [nchain x N struct]: fields are "x" for the parameter vector
%      and "f" for the value of the parameter vector
%      - **bestf** [numeric]: best function value
%      - **bestx** [vector]: best parameter vector
%      - **lambda** [vector]: tempering levels
%      - **best** [nchain x 1]: vector of best individual in each chain
%      - **m** [vector]: mean of the parameter draws
%      - **SIG** [matrix]: covariance of the parameter draws
%      - **m_algo** [vector]: mean with particular updating
%      - **SIG_algo** [matrix]: covariance with particular updating
%      - **funevals** [integer]: function evaluations
%      - **log_mdd** [numeric]: log marginal data density
%      - **ess** [numeric]: effective sample size
%      - **individual_cores_stats** [struct]: stats on the optimization
%
% Note:
%
%    - It is assumed that logf is a function to minimize
%
%    - the update of c in one bin might not be adequate for the other bins. As
%    a result, the last bin may be detrimental for the first bin. This is
%    calls for bin-specific log_c: IMPLEMENTED!!!
%
%    - Possible variations for the stud
%      - best in the bin (Current implementation)
%      - random
%      - simple average : what if it has unacceptable density?
%      - weigthed average through a merit system : what if it has unacceptable density?
%
% Example:
%
%    See also:	CONSTANT_BVAR_SAMPLER, MH_SAMPLER

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
typical_nshocks=20;
typical_sample=200;
defaults={ % arg_names -- defaults -- checks -- error_msg
    'M_over_N',5/100,@(x)num_fin(x) && x>0 && x<1,'M_over_N (percentage per bin) must be in (0,1)'
    'c',1,@(x)num_fin(x) && x>0,'c (tuning parameter) should be a positive scalar'
    'c_range',[sqrt(eps),100],@(x)numel(x)==2 && all(num_fin(x) & x>0) && x(2)>=x(1),'c_range variation range for tuning parameter must be a two element vector'
    'fixed_scaling',false,@(x)isscalar(x) && islogical(x),'fixed_scaling (fixed tuning parameter) should be a logical scalar'
    'use_true_moments',false,@(x)isscalar(x) && islogical(x),'use_true_moments should be a logical scalar'
    'N',5000,@(x)num_fin(x) && x>0,'N(# simulations per core) should be a strictly positive integer'
    'H',40,@(x)isempty(x)||(num_fin(x) && x>0),'H(# tempering stages) should be a strictly positive integer or empty'
    'p',[],@(x)isempty(x)||(num_fin(x) && x>=0 && x<1),'p (prob of drawing from prev. distrib.) must be in [0,1) when not empty. If empty, it will be set to 0.1*ps'
    'lambda_1',1/(10*typical_nshocks*typical_sample),@(x)num_fin(x) && x>0 && x<=1,'lambda_1 (first level of tempering) should be in (0,1]'
    'geometric_lambda',true,@(x)isscalar(x) && islogical(x),'geometric_lambda should be a logical scalar'
    'ps',.9,@(x)num_fin(x) && x>0 && x<=1,'ps (prob of storing a metropolis draw) must be in (0,1]'
    'alpha',[.25,.45],@(x)((numel(x)==2 && x(2)>=x(1))||numel(x)==1) && all(x>0) && all(x<1),'target_range should be a 2-element vector s.t. x(2)>=x(1) and x(i) in (0,1)'
    'rwm_exp',0.6,@(x)num_fin(x) && x>0.5 && x<1,'rwm_exp (Exponent of random-walk adaptation step size) must be in (1/2,1)'
    'ess_min',0.1,@(x)num_fin(x) && x>0 && x<1, 'ess_min (constant for controlling the evolution of lambda) should be in (0,1)'
    'p_mutant',0,@(x)num_fin(x) && x>=0 && x<=1,'p_mutant (prob of drawing a mutant) must be in [0,1]'
    'penalty',1e+8,@(x)isempty(x)||num_fin(x) && x>0,'penalty (worst value possible in absolute value) must be empty or a finite positive number'
    'extra_runs',5,@(x)num_fin_int(x) && x>0,' extra_runs (extra runs after reaching tempering of 1) must be a positive integer'
    'focus_fire_after',.5,@(x)num_fin(x) && x>=0 && x<=1,'focus_fire_after (focus on the best after) must be in [0,1]'
    };

if nargin==0
    Results=cell2struct(defaults(:,2),defaults(:,1),1);
    return
end

if nargin<6
    SIG=[];
    if nargin<5
        mu=[];
        if nargin<4
            options=struct();
        end
    end
end

% number of parameters
[d,ncols]=size(lb);
if ncols~=1
    error('number of columns of lb should be 1')
end

if ~(all(isfinite(lb)) && all(isfinite(ub)))
    error('lb and ub shoud be finite')
end

if d==0
    error('lb cannot have less than 1 row')
end

if ~isequal(size(ub),[d,ncols])
    error('size ub does not match size lb')
end

options=parse_arguments(defaults,options);

if ischar(logf)
    logf=str2func(logf);
end

if ~isa(logf,'function_handle')
    error('logf should be a function handle or a string')
end

% worst value that the function can assume
%------------------------------------------
if isempty(options.penalty)
    options.penalty=1e+6;
end
penalty=options.penalty;

% number of striations
%----------------------
% choose L so that the there is an equal number of draws in each striation
M_over_N=options.M_over_N;
striations=0:M_over_N:1;
M=numel(striations)-1;

% tuning parameter for the covariance of the metropolis
%-------------------------------------------------------
log_c=log(options.c)*ones(1,M);

% adaptation
%-------------
rwm_exp=options.rwm_exp;
alpha=options.alpha;
fixed_scaling=options.fixed_scaling;
c_range=options.c_range;
use_true_moments=options.use_true_moments;

% number of simulations per core
%-----------------------
N=options.N;
number_of_draws_per_striation=round(N/M);

if isempty(mu)
    mu=.5*(ub+lb);
end

if isempty(SIG)
    SIG=1e-4*eye(d);%SIG=covar;
end

mu_algo=mu;
SIG_algo=SIG;

% draw initial distribution:
%----------------------------
[genpop,funevals,log_I_old,~,NumWorkers]=utils.mcmc.initial_draws(logf,lb,ub,N,penalty,mu);
% funevals is a vector so squash it right here right now
%-------------------------------------------------------
funevals=sum(funevals);
bestf=genpop(1).f;

% focus fire after
%------------------
focus_fire_after=options.focus_fire_after;

is_time_to_exploit=false;

% Tempering levels
%------------------
geometric_lambda=options.geometric_lambda;

% number of stages
%------------------
H=options.H;
fixed_horizon=~isempty(H);

lambda_vector=tempering_steps(options.lambda_1);

% probability of storing a metropolis draw
%-------------------------------------------
ps=options.ps;

% probability of drawing from the distribution of the previous stage.
%----------------------------------------------------------------------
% Instead of just drawing from the previous stage, one can devise a
% bee-type of algorithm within the striation and still use the metropolis
% to accept or reject
if isempty(options.p)
    options.p=0.1*ps;
end
total_runs=H+options.extra_runs;
log_I=cell(1,H);

istage=0;
done=false;
accept_ratio=[];
n_vectors=[];
nbatch=0;
last_generation=genpop;
% create a file that will abort estimation if too long
utils.optim.manual_stopping();
I_am_done=false;
while ~done
    istage=istage+1;
    
    is_time_to_exploit=fixed_horizon && istage>=focus_fire_after*H;
    
    lambda_old=lambda_vector(istage);
    
    lambda=lambda_vector(istage+1);
    
    [genpop,log_I{istage},Fnewpop]=do_one_stage(logf,genpop,log_I_old);
    
    % compute log marginal data density 
    %------------------------------------    
    log_I_old=log_I{istage};
    
    if fixed_horizon
        done=istage==total_runs;
        % the lambda is already computed
    else
        done=lambda_converged(lambda);
        if ~done
            % update the grand lambda at least for the computation of overall
            % moments even if the individual cores use specific moments
            %-----------------------------------------------------------------
            lambda_vector(istage+2)=find_next_lambda(lambda,Fnewpop,...
                options.minimization,options.ess_min);
        end
    end
    done=done||I_am_done;
end

% remove the common options in each suboption
%---------------------------------------------
fields=fieldnames(options);
trimmed_options=rmfield(options,fields);

Results=struct('pop',genpop,...
    'bestf',genpop(1).f,...
    'bestx',genpop(1).x,...
    'lambda',lambda_vector,...
    'm',mu,...
    'SIG',SIG,...
    'm_algo',mu_algo,...
    'SIG_algo',SIG_algo,...
    'funevals',funevals,...
    'log_mdd',cell2mat(log_I),...
    'ess',options.ess,...
    'individual_cores_stats',trimmed_options);

    function lambda_vector=tempering_steps(lambda_1)
        if fixed_horizon
            if geometric_lambda
                b=(1/lambda_1)^(1/(H-1));
                lambda_vector=lambda_1*b.^(0:H-1);
            else
                lambda_vector=linspace(lambda_1,1,H);
            end
            lambda_vector=[0,lambda_vector];
            if abs(lambda_vector(end)-1)>1e-9
                error('error in the computation of the tempering steps')
            end
            lambda_vector(end)=1;
            lambda_vector=[lambda_vector,ones(1,options.extra_runs)];
        else
            lambda_vector=[0,lambda_1];
        end
    end

    function [pop,log_I,Fnewpop]=do_one_stage(logf,pop0,log_I_old)
        p=options.p;
        ps=options.ps;
        probability_of_not_storing_metropolis=1-ps;
        p_mutant=options.p_mutant;
        % handle to a key function
        %--------------------------
        dsmh_draw_hdle=@dsmh_draw;
        
        if is_time_to_exploit
            target_pop=last_generation;
        else
            target_pop=pop0;
        end
        
        % it is assumed the population is already sorted
        %------------------------------------------------
        F_pop=[target_pop.f];
        
        % separate in striations: here we use the actual number of
        % individuals instead of N
        %------------------------------------------------------------
        separators=round(striations*numel(F_pop));
        ndraws_per_striation=separators(2:end)-separators(1:end-1);
        target_pop=mat2cell(target_pop,1,ndraws_per_striation);
        d=size(lb,1);
        pop=cell(1,M);
        
        tmp_all=tic;
        loop_is_complete=false(1,M);
        nbatch=nbatch+1;
        metrop_draws=zeros(1,M);
        metrop_accept=zeros(1,M);
        parfor (istri=1:M,NumWorkers)
%             if I_should_stop
%                 % basically check convergence
%                 continue
%             end
            tmp_i=tic;
            if use_true_moments
                sqrt_cSIG=chol(exp(log_c(istri))*SIG,'lower');
            else
                sqrt_cSIG=chol(exp(log_c(istri))*SIG_algo,'lower');
            end
            pointer=istri;
            theta_star=target_pop{pointer}(1);
            pop{istri}=theta_star(1,ones(1,number_of_draws_per_striation));
            for idraw=1:number_of_draws_per_striation
                [theta_star,metrop_accept(istri),metrop_draws(istri)]=dsmh_draw_hdle(...
                    theta_star,metrop_accept(istri),metrop_draws(istri),...
                    sqrt_cSIG,pointer);
                pop{istri}(idraw)=theta_star;
            end
            
            accept_ratio=nan;
            if metrop_draws(istri)
                % note we use the same adaptation exponent
                % (options.rwm_exp) both for the moments and for the scale
                % of the covariance matrix.
                accept_ratio=metrop_accept(istri)/metrop_draws(istri);
                log_c(istri)=utils.mcmc.update_scaling(log_c(istri),...
                    accept_ratio,alpha,fixed_scaling,...
                    nbatch,rwm_exp,[],c_range);
            end
            % N.B: Two opposing forces: A flat tempering makes it easy to
            % accept. When we accept too much we make the covariance matrix
            % bigger and so we accept vectors far away. By the time the
            % tempering gets to more normal levels, the covariance matrix
            % is perhaps too big to decrease quickly enough. This may call
            % for a geometric lambda
            
            % display partial results
            %------------------------
            pop_i=utils.mcmc.sort_population(pop{istri});
            I_am_done=display_progress(istage,total_runs,istri,M,lambda,funevals,...
                pop_i(1).f,accept_ratio,[],[],toc(tmp_i),is_time_to_exploit);
            loop_is_complete(istri)=true;
        end
        options.metrop_draws=sum(metrop_draws);
        options.metrop_accept=sum(metrop_accept);
        % resort the new database
        %-------------------------
        % add index in case the procedure is stopped before the loop was
        % completed.
        pop=utils.mcmc.sort_population(cell2mat(pop(loop_is_complete)));
        
        % update moments to be used in the metropolis step and the
        % theoretical moments
        %-------------------------------------------------------------
        [mu_algo,SIG_algo]=utils.mcmc.update_moments(mu_algo,SIG_algo,...
            [pop.x],nbatch,rwm_exp);
        
        [mu,SIG,n_vectors]=utils.moments.recursive(mu,SIG,[pop.x],n_vectors);
        % correct covariance if necessary
        %--------------------------------
        SIG=utils.cov.project(SIG);
            
        options.accept_ratio=options.metrop_accept/options.metrop_draws;
        
        last_generation=pop;
        
        Fnewpop=[pop.f];
        % calculate some statistics on the new population
        %-------------------------------------------------
        [w,wtilde,options.ess]=importance_weights(Fnewpop,lambda,lambda_old);
        
        log_I=log_I_old+log(mean(wtilde(:)));
        
        % merge populations
        %------------------
        pop=[pop,pop0];
        
        % resort the database
        %---------------------
        pop=utils.mcmc.sort_population(pop);
        
        bestf=pop(1).f;
        
        time_it_took=toc(tmp_all);
        I_am_done=display_progress(istage,total_runs,[],[],lambda,funevals,...
            bestf,options.accept_ratio,options.ess,numel(w),time_it_took,...
            is_time_to_exploit);
        
        function [theta,metrop_accept__,metrop_draws__]=dsmh_draw(theta_star,...
                metrop_accept__,metrop_draws__,sqrt_cSIG,pointer)
            
            is_metrop_draw=rand>p;
            if is_metrop_draw
                % gaussian distribution
                %--------------------------
                [theta,lgratio]=gaussian_draw();
            else
                % previous level within the same striation
                %-------------------------------------------
                [theta,lgratio]=previous_level_draw();
            end
            
            [theta_star]=dsmh_selection(theta_star,theta);
            
            % dynamic thinning
            %------------------
            while is_metrop_draw && rand<probability_of_not_storing_metropolis
                [theta,lgratio]=gaussian_draw();
                [theta_star]=dsmh_selection(theta_star,theta);
            end
            
            % update final
            %--------------
            theta=theta_star;
            
            function [theta,accepted]=dsmh_selection(theta0,theta)
                if utils.mcmc.is_violation(theta.f,penalty)
                    theta=theta0;
                    accepted=false;
                else
                    lfratio=lambda*(theta.f-theta0.f);
                    lfg_ratio=lfratio+lgratio;
                    
                    % minimization: change sign
                    %--------------------------
                    lfg_ratio=-lfg_ratio;
                    
                    % Metropolis criterion
                    %-----------------------
                    crit=min(0,lfg_ratio);
                    accepted=log(rand)<crit;
                    if accepted
                        if is_metrop_draw
                            metrop_accept__=metrop_accept__+1;
                        end
                    else
                        theta=theta0;
                    end
                end
            end
            
            function [theta,lgratio]=previous_level_draw()
                % find the right striation for the newcomer
                target_stri=pointer;
                nail=randi(numel(target_pop{target_stri}));
                
                if rand<p_mutant
                    [theta]=draw_mutant();
                else
                    theta=target_pop{target_stri}(nail);
                end
                lgratio=lambda_old*(theta_star.f-theta.f);
                
                function [theta]=draw_mutant()
                    xd=generate_mutant([target_pop{target_stri}.x],nail,lb,ub);
                    theta.f=logf(xd);
                    theta.x=xd;
                    funevals=funevals+1;
                end
            end
            
            function [theta,lgratio]=gaussian_draw()
                xd=theta_star.x+sqrt_cSIG*randn(d,1);
                bad=xd<lb; xd(bad)=lb(bad);
                bad=xd>ub; xd(bad)=ub(bad);
                theta.f=logf(xd);
                theta.x=xd;
                funevals=funevals+1;
                lgratio=0;
                metrop_draws__=metrop_draws__+1;
            end
        end
    end

end

function mutant=generate_mutant(xx,ii,lb,ub)
[d,nx]=size(xx);
% generate a solution for index ii
% a randomly chosen solution different from ii is used
% for producing a mutant
while 1
    donor_id=min(nx,fix(rand*nx)+1);
    if donor_id~=ii
        break
    end
end
mutant=xx(:,ii);

% pick the parameter to change in the new solution
%--------------------------------------------------
change=min(fix(rand*d)+1,d);
nc=numel(change);
mutant(change)=mutant(change)+(mutant(change)-xx(change,donor_id))*2.*(rand(nc,1)-.5);

mutant(change)=utils.optim.recenter(mutant(change),lb(change),ub(change));

end

function flag=lambda_converged(lambda)
flag=1-lambda<sqrt(eps);
end

function I_am_done=display_progress(istage,nstages,bin,nbins,lambda,funevals,bestf,...
    accept_ratio,ess,N,nsecs,is_time_to_exploit)
if nargin<10
    nsecs=nan;
end
if is_time_to_exploit
    typeof='climbing';
else
    typeof='exploring';
end
if isempty(bin)
    string=['stage %0.0f/%0.0f(%s), lambda %8.4f, f-count %4.0f, best f(x) %8.4f'...
        ', accept. rate %4.2f, ESS(N) %4.4f(%0.0f), this iteration %0.4f sec'];
    data={istage,nstages,typeof,lambda,funevals,bestf,100*accept_ratio,ess,N,nsecs};
else
    string=['stage %0.0f/%0.0f(%s), bin %0.0f/%0.0f, lambda %8.4f'...
        ', f-count %4.0f, best f(x) %8.4f, accept. rate %4.2f, this iteration %0.4f sec'];
    data={istage,nstages,typeof,bin,nbins,lambda,funevals,bestf,100*accept_ratio,nsecs};
end

fprintf([string,'\n\n'],data{:});

I_am_done=how_to_be_done();

    function I_am_done=how_to_be_done()
        I_am_done=false;
        ms=utils.optim.manual_stopping(1);
        if ms~=0
            if ms==2
                keyboard
            else
                I_am_done=true;
                disp('manually aborted (without extreme prejudice)')
            end
        end
    end
end

function [lamb,fval,exitflag,output,jacobian]=find_next_lambda(lambda_old,F,ess_min,solver)
if nargin<5
    solver='lsqnonlin';
    if nargin<4
        ess_min=0.1;
    end
end
N=numel(F);
N_ess_min=N*ess_min;

% start in the middle
%--------------------
x0=.5*(lambda_old+1);
options=optimset('Display','none');
switch solver
    case 'lsqnonlin'
        [lamb,resnorm,~,exitflag,output,~,jacobian] = lsqnonlin(...
            @objective,x0,lambda_old,1,options);
        fval=resnorm;
    case 'fsolve'
        [lamb,fval,exitflag,output,jacobian] =fsolve(@objective,x0,options);
    case 'fzero'
        [lamb,fval,exitflag,output] = fzero(@objective,x0,options);
    otherwise
        error(['unknown solver ',solver])
end
if ~(abs(fval)<sqrt(eps)) && (1-lamb>sqrt(eps))
    % if lambda==1, then fval could fail but if lamb==1 we have reached the
    % final level of tempering and this is ok
    error('finding lambda failed')
end
    function discrep=objective(lambda)
        % implement bounds for fsolve and fzero
        %---------------------------------------
        if lambda<lambda_old||lambda>1
            discrep=1e+8;
            if lambda>1
                discrep=-discrep;
            end
        else
            [~,~,ess]=importance_weights(F,lambda,lambda_old);
            discrep=ess-N_ess_min;
        end
    end
end

function [w,wtilde,ess]=importance_weights(F,lambda,lambda_old)
% the weights may collapse if the lambda differential is so wide that in
% case F has large components, the exponential goes to infinity. This
% probably happens because the initial distribution is not a probability
% density...
%-----------------------------------------------------------------------

% minimization: change sign
%--------------------------
F=-F;

wtilde=(lambda-lambda_old)*F;
% log_big=max(wtilde);
wtilde=exp(wtilde);% wtilde=exp(wtilde-log_big);
w=wtilde/sum(wtilde);
if nargout>2
    % effective sample size
    %------------------------
    ess=effective_sample_size(w);
end
if any(isnan(wtilde))||any(isinf(wtilde))
    save problematic_importance_weights
end
end

function ess=effective_sample_size(w)
ess=1/sum(w(:).^2);
end