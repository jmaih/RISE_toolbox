function [Results]=mh_sampler(logf,lb,ub,options,mu,SIG)
% MH_SAMPLER -- Metropolis Hastings sampler
%
% Syntax
% -------
% ::
%
%   [Results]=MH_SAMPLER(logf,lb,ub)
%
%   [Results]=MH_SAMPLER(logf,lb,ub,options)
%
%   [Results]=MH_SAMPLER(logf,lb,ub,options,mu)
%
%   [Results]=MH_SAMPLER(logf,lb,ub,options,mu,SIG)
%
% Inputs
% -------
%
% - **logf** [char|function_handle]: Objective function to MINIMIZE!!!
%
% - **lb** [d x 1 vector]: lower bound of the paramters
%
% - **ub** [d x 1 vector]: upper bound of the paramters
%
% - **options** [struct]:
%   - **alpha** [scalar|2-element|{[.25,.45]}]: target acceptance rate
%   - **burnin** [integer|{0}]: number of burn-in initial simulations
%   - **N** [integer|{20000}]: number of simulations
%   - **verbose** [integer|{100}]: displays progress for every multiple of
%   "verbose"
%   - **c** [scalar|{0.25}]: initial scale for the covariance matrix
%   - **c_range** [vector|{}]: range of variation of c
%   - **thin** [integer|{1}]: number of thinning simulations. 1 means we
%   keep every draw, 2 means we keep every second, 3 every third, etc.
%   - **retune_cov_every** [integer|{100}]: frequence for the retuning of
%   the scale parameter
%   - **penalty** [numeric|{1e+8}]: worst possible function value
%   - **nchain** [integer|{1}]: number of parallel chains
%   - **rwm_exp** [numeric|{0.6}]: tuning hyper-parameter for scale and
%   covariance matrix
%   - **fixed_scaling** [true|{false}]: if true, the scaling (c) of the
%   covariance matrix is kept constant
%   - **use_true_moments** [true|{false}]: if true, the updated exact
%   covariance matrix of the draws is used at each step. If false, a
%   different update of the covariance matrix is used.
%   - **logproppdf** [function_handle|{[]}]: used when the proposal is not
%   symmetric
%   - **MaxTime** [numeric|{inf}]: maximum simulation time.
%   - **adapt_covariance** [true|{false}]: If true, the covariance matrix
%   is updated with the sampled draws.
%
% - **mu** [d x 1 vector]: initial condition for the sampler
%
% - **SIG** [d x d matrix]: initial covariance matrix
%
% Outputs
% --------
%
% - **Results** [struct]:
%   - **pop** [nchain x N struct]: fields are "x" for the parameter vector
%   and "f" for the value of the parameter vector
%   - **bestf** [numeric]: best function value
%   - **bestx** [vector]: best parameter vector
%   - **best** [nchain x 1]: vector of best individual in each chain
%   - **m** [vector]: mean of the parameter draws
%   - **SIG** [matrix]: covariance of the parameter draws
%   - **m_algo** [vector]: mean with particular updating
%   - **SIG_algo** [matrix]: covariance with particular updating
%   - **funevals** [1 x nchain vector]: function evaluations
%   - **stats** [struct]: stats on the optimization
%
% More About
% ------------
%
% - It is assumed that logf is a function to minimize
%
% Examples
% ---------
%
% See also:	CONSTANT_BVAR_SAMPLER, RRF_SAMPLER

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
defaults={ % arg_names -- defaults -- checks -- error_msg
    'alpha',[.25,.45],@(x)((numel(x)==2 && x(2)>=x(1))||numel(x)==1) && all(x>0) && all(x<1),'target_range should be a 2-element vector s.t. x(2)>=x(1) and x(i) in (0,1)'
    'burnin',0,@(x)num_fin_int(x),'burnin should be an integer in [0,inf)'
    'N',20000,@(x)num_fin_int(x) && x>0 ,'N should be a strictly positive integer'
    'verbose',100,@(x)num_fin_int(x) && x>0 ,'verbose should be a strictly positive integer'
    'c',0.25,@(x)num_fin(x) && x>0 ,'c (tuning parameter) should be a positive scalar'
    'c_range',[sqrt(eps),100],@(x)numel(x)==2 && all(num_fin(x) & x>0) && x(2)>=x(1),'c_range variation range for tuning parameter must be a two element vector'
    'thin',1,@(x)num_fin_int(x) && x>=1,'thin must be >=1'
    'retune_cov_every',100,@(x)num_fin_int && x>0, 'retune_cov_every should be a positive integer'
    'penalty',1e+8,@(x)isempty(x)||num_fin(x) && x>0,'penalty (worst value possible in absolute value) must be empty or a finite positive number'
    'nchain',1,@(x)num_fin_int(x) && x>0,'nchain(# chains) should be a strictly positive integer'
    'rwm_exp',0.6,@(x)num_fin(x) && x>0.5 && x<1,'rwm_exp (Exponent of random-walk adaptation step size) must be in (1/2,1)'
    'fixed_scaling',false,@(x)isscalar(x) && islogical(x),'fixed_scaling (fixed tuning parameter) should be a logical scalar'
    'use_true_moments',false,@(x)isscalar(x) && islogical(x),'use_true_moments should be a logical scalar'
    'logproppdf',[],@(x)isa(x,'function_handle'),'logproppdf must be a function handle'
    'MaxTime',inf,@(x)num_fin(x) && x>0,'MaxTime must be a positive scalar'
    'adapt_covariance',false,@(x)isscalar(x) && islogical(x),'adapt_covariance should be a logical scalar'
%     'delay_rejection',true,@(x)isscalar(x) && islogical(x),'delay_rejection should be a logical scalar'
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

logproppdf=options.logproppdf;

symmetric=isempty(logproppdf);

if ~isa(logf,'function_handle')
    error('logf should be a function handle or a string')
end

% target acceptance range
%------------------------
alpha=options.alpha;

% more options
%-------------
fixed_scaling=options.fixed_scaling;
rwm_exp=options.rwm_exp;
c_range=options.c_range;
use_true_moments=options.use_true_moments;

% worst value that the function can assume
%------------------------------------------
penalty=options.penalty;

% burn-in
%--------
burnin=options.burnin;

% thining
%--------
thin=options.thin;

% number of chains
%------------------
nchain=options.nchain;

% tuning parameter for the covariance of the metropolis
%-------------------------------------------------------
log_c=log(options.c)*ones(1,nchain);

% number of simulations per chain
%----------------------------------
N=options.N;

if isempty(mu)
    mu=.5*(ub+lb);
end

% covariance adaptations
%-------------------------
adapt_covariance=options.adapt_covariance;

if isempty(SIG)
    SIG=1e-4*eye(d);%SIG=covar;
    adapt_covariance=true;
end
mu=mu(:,ones(1,nchain));
SIG=SIG(:,:,ones(1,nchain));
mu_algo=mu;
SIG_algo=SIG;

% draw initial distribution:
%----------------------------
[stud,funevals,~,~,NumWorkers]=utils.mcmc.initial_draws(logf,lb,ub,nchain,penalty,mu);
% vectorize in case of many parallel chains
stud=stud(:);

% pre-allocate
%--------------
pop=stud(:,ones(1,N));

retune_cov_every=options.retune_cov_every;

d=numel(lb);

obj=struct();
obj.funcCount=sum(funevals);
obj.iterations=0;
obj.start_time=clock;
obj.MaxTime=options.MaxTime;
obj.verbose=options.verbose;
obj.MaxFunEvals=inf;
obj.best_fval=[];
obj.MaxNodes=1;
obj.optimizer=mfilename;
obj.number_of_parameters=size(mu,1);
obj.accept_ratio=[];

sqrt_cSIG=[];
nbatch=1;
accept=zeros(1,nchain);
big_accept=0;
scale_times_sqrt_covariance_updating()

idraw=-burnin;
total_draws=N*thin+burnin;
obj.MaxIter=total_draws;
U = log(rand(total_draws,nchain));
q_x0_given_y = 0;
q_y_given_x0 = q_x0_given_y;

utils.optim.manual_stopping;

stopflag=utils.optim.check_convergence(obj);

best=stud;

obj.best_fval=[best.f];

accept_ratio=[];

% wtbh=waitbar(0,'please wait...','Name','MH sampling');

while isempty(stopflag)
    idraw=idraw+1;
    obj.iterations=idraw+burnin;
    % sample from the proposal distribution
    %---------------------------------------
    y=proprnd();
    if ~symmetric
        q_x0_given_y = logproppdf(stud,y);
        q_y_given_x0 = logproppdf(y,stud);
    end
    % this is a generic formula.
    rho = ([y.f]+q_x0_given_y)-([stud.f]+q_y_given_x0);
    % minimization: change sign
    %--------------------------
    rho=-rho;
    
    % Accept or reject the proposal
    %-------------------------------
    acc = min(rho,0)>=U(obj.iterations,:);
    
    stud(acc,:) = y(acc,:); % preserves x's shape.
    
    accept = accept+acc;
    
    % save down
    %----------
    if idraw>0 && mod(idraw,thin)==0 % burnin and thin
        
        the_iter=idraw/thin;
        
        pop(:,the_iter) = stud;
        
    end
    
    % update the best
    %-----------------
    tmpf=[stud.f];
    
    good=tmpf<obj.best_fval;
    
    best(good)=stud(good);
    
    obj.best_fval=[best.f];
    
    % update the moments
    %-------------------
    moments_updating()
    
    % update the cCs
    %---------------
    if mod(obj.iterations,retune_cov_every)==0
        nbatch=nbatch+1;
        scale_times_sqrt_covariance_updating()
    end
    
    % display progress
    %------------------
    obj.accept_ratio=accept_ratio;
    utils.optim.display_progress(obj)
    
    % check convergence
    %------------------
    stopflag=utils.optim.check_convergence(obj);
    
%     waitbar_updating()
end

% delete(wtbh)

do_results()

    function do_results()
        pop=pop(:,1:the_iter);
        bestOfTheBest=utils.mcmc.sort_population(best);
        obj.end_time=clock;
        Results=struct('pop',pop,...
            'bestf',bestOfTheBest(1).f,...
            'bestx',bestOfTheBest(1).x,...
            'best',best,...
            'm',mu,...
            'SIG',SIG,...
            'm_algo',mu_algo,...
            'SIG_algo',SIG_algo,...
            'funevals',funevals,...
            'stats',obj);
    end

    function scale_times_sqrt_covariance_updating()
        accept_ratio=accept/retune_cov_every;
        sqrt_cSIG=SIG;
        iterations=obj.iterations;
        parfor (ichain=1:nchain,NumWorkers)
            if iterations>0
                log_c(ichain)=utils.mcmc.update_scaling(log_c(ichain),...
                    accept_ratio(ichain),alpha,fixed_scaling,nbatch,...
                    rwm_exp,[],c_range);
            end
            if use_true_moments
                sqrt_cSIG(:,:,ichain)=chol(...
                    exp(log_c(ichain))*SIG(:,:,ichain),...
                    'lower');
            else
                sqrt_cSIG(:,:,ichain)=chol(...
                    exp(log_c(ichain))*SIG_algo(:,:,ichain),...
                    'lower');
            end
        end
        big_accept=big_accept+accept;
        accept=0*accept;
    end

    function moments_updating()
        if adapt_covariance
            for ichain=1:nchain
                [mu(:,ichain),SIG(:,:,ichain)]=utils.moments.recursive(...
                    mu(:,ichain),SIG(:,:,ichain),stud(ichain).x,obj.iterations);
                [mu_algo(:,ichain),SIG_algo(:,:,ichain)]=utils.mcmc.update_moments(...
                    mu_algo(:,ichain),SIG_algo(:,:,ichain),stud(ichain).x,obj.iterations,...
                    options.rwm_exp);
            end
        end
    end

    function v1=proprnd()
        v1=stud;
        parfor (ichain=1:nchain,NumWorkers)
            xd=stud(ichain).x+sqrt_cSIG(:,:,ichain)*randn(d,1);
            bad=xd<lb; xd(bad)=lb(bad);
            bad=xd>ub; xd(bad)=ub(bad);
            v1(ichain).x=xd;
            v1(ichain).f=logf(xd); %#ok<PFBNS>
        end
        funevals=funevals+nchain;
        obj.funcCount=sum(funevals);
    end

%     function waitbar_updating()
%         x=obj.iterations/total_draws;
%         waitbar(x,wtbh,...
%             {
%             sprintf('bestf %s',num2str(obj.best_fval))
%             sprintf('acceptance rate %s',num2str(100*accept_ratio))
%             }...
%             )
%     end

end
