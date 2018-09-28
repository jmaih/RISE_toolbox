function Results = adaptive_parallel_tempering(logf,lb,ub,options,mu,SIG)
% Also known as:
% - Replica exchange Monte Carlo
% - Metropolis coupled Markov chain Monte Carlo
%
% Out:
%   X      -- Simulated variables in 3D array (dim*levels*N)
%   m      -- The final proposal means (vectors of length dim)
%   Rchol      -- The final proposal covariances (dim*dim matrices)
%   log_c  -- The final scaling factors (vector of length H)
%   lambda   -- The final inverse temperatures (vector of length H)
%   stats  -- Acceptance rate statistics etc.

num=@(x)isnumeric(x) && isscalar(x) && isreal(x);
num_fin=@(x)num(x) && isfinite(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
% typical_nshocks=20;
% typical_sample=200;
defaults={ % arg_names -- defaults -- checks -- error_msg
    'N',5000,@(x)num_fin_int(x) && x>0,'N(# loops/cycles) must be a positive integer'
    'H',12,@(x)num_fin(x) && x>0,'H(# tempering stages) should be a strictly positive integer'
    'minimization',true,@(x)isscalar(x) && islogical(x),'minimization must be a logical scalar'
    'c',1,@(x)num_fin(x) && x>0,'c (cov tuning parameter) should be a real and positive scalar'
    'alpha',.234,@(x)num_fin(x) && x>0 && x<=1,'alpha (target acceptance rate) must be in (0,1]'
    'alpha_swap',.234,@(x)num_fin(x) && x>0 && x<=1,'alpha_swap (target acceptance rate for swaps) must be in (0,1]'
    'nthin',1,@(x)num_fin_int(x) && x>0,'nthin(# thinning draws) must be a positive integer'
    'verbose_iter',10,@(x)num_fin_int(x),'verbose_iter (show every verbose_iter iteration) must be a positive integer'
    'rwm_fixed_p',0,@(x)num_fin(x) && x>=0 && x<=1,'rwm_fixed_p (prob of drawing from a fixed initial proposal) must be in [0,1]'
    'rwm_exp',0.6,@(x)num_fin(x) && x>=0,'rwm_exp (Exponent of random-walk adaptation step size) must be a positive scalar'
    'sw_exp_rm',0.6,@(x)num_fin(x) && x>=0,'sw_exp_rm (Exp. of temperature adaptation step size) must be a positive scalar'
    'ram_adapt',false,@(x)isscalar(x) && islogical(x),'ram_adapt (Use the robust AM adaptation) must be a logical scalar'
    'separate_shape_adaptation',true,@(x)isscalar(x) && islogical(x),'separate_shape_adaptation (separate covariance for each temperature) must be a logical scalar'
    'MaxFunEvals',inf,@(x)isempty(x)||num(x) && floor(x)==ceil(x) && x>0,'MaxFunEvals (max # func evals) must be a positive integer'
    'MaxTime',inf,@(x)isempty(x)||num(x) && x>0,'MaxTime (max time) must be a positive scalar when not empty'
    'fixed_temperatures',false,@(x)isscalar(x) && islogical(x),'fixed_temperatures (fixed temp) must be a logical scalar'
    'penalty',1e+8,@(x)isempty(x)||num_fin(x) && x>0,'penalty (worst value possible in absolute value) must be empty or a finite positive number'
    'debug',false,@(x)isscalar(x) && islogical(x),'debug must be a logical scalar'
    'geometric_tempering',true,@(x)isscalar(x) && islogical(x),'geometric_tempering must be a logical scalar'
    'rng',[],@(x)isscalar(x),'rng must be empty or a scalar'
    };
    %     'nu',5,@(x)num_fin_int(x) && x>0,' nu must be a positive integer'
    %     'M',20,@(x)num_fin_int(x) && x>=1,'M (# striations) must be a positive integer'
    %     'recombine_all',false,@(x)isscalar(x) && islogical(x),'recombine_all should be a logical scalar'
    %     'NG',5000,@(x)num_fin(x) && x>0,'NG(# simulations) should be a strictly positive integer'
    %     'G',1,@(x)num_fin(x) && x>0,'G(# cores) should be a strictly positive integer'
    %     'p',[],@(x)isempty(x)||(num_fin(x) && x>=0 && x<1),'p (prob of drawing from prev. distrib.) must be in [0,1) when not empty. If empty, it will be set to 0.1*ps'
    %     'lambda_1',1/(10*typical_nshocks*typical_sample),@(x)num_fin(x) && x>0 && x<=1,'lambda_1 (first level of tempering) should be in (0,1]'
    %     'geometric_lambda',true,@(x)isscalar(x) && islogical(x),'geometric_lambda should be a logical scalar'
    %     'ps',.9,@(x)num_fin(x) && x>0 && x<=1,'ps (prob of storing a metropolis draw) must be in (0,1]'

% % degrees of freedom of multivariate t-student distribution
% %-----------------------------------------------------------
% nu=options.nu;

% % number of striations
% %----------------------
% M=options.M;

% ess_iw_min=0.1;

% % recombination of particiles/bees/vectors
% %-------------------------------------------
% % if recombine_all==false, we leave each striation unchanged... the only
% % communication with the other striations is through the common covariance
% % matrix
% recombine_all=options.recombine_all;
% % number of cores
% %-----------------
% % default should be 1
% G=options.G;

% % probability of storing a metropolis draw
% %-------------------------------------------
% ps=options.ps;

% % probability of drawing from the distribution of the previous stage.
% %----------------------------------------------------------------------
% % Instead of just drawing from the previous stage, one can devise a
% % bee-type of algorithm within the striation and still use the metropolis
% % to accept or reject
% p=options.p;
% if isempty(p)
%     p=0.1*ps;
% end
% % Tempering levels
% %------------------
% geometric_lambda=options.geometric_lambda;
%
% lambda_vector=tempering_steps(options.lambda_1);

if nargin==0
    Results=cell2struct(defaults(:,2),defaults(:,1),1);
    return
end

if nargin<6
    SIG=[];
    if nargin<5
        mu=[];
    end
end

options=utils.miscellaneous.parse_arguments(defaults,options);

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

if ischar(logf)
    logf=str2func(logf);
end

if ~isa(logf,'function_handle')
    error('logf should be a function handle or a string')
end

if ~isempty(options.rng)
    rng(options.rng)
end

% worst value that the function can assume
%------------------------------------------
penalty=options.penalty;

% debuging
%-----------
debug=options.debug;

geometric_tempering=options.geometric_tempering;

% fixed or time-varying temperatures
%------------------------------------
fixed_temperatures=options.fixed_temperatures;

% option for whether the optimization is a maximization or a minimization
% ------------------------------------------------------------------------
minimization=options.minimization;

% number of stages
%------------------
H=options.H;

% tuning parameter for the covariance of the metropolis
%-------------------------------------------------------
log_c = log(options.c);
if numel(log_c) == 1
    log_c=log_c*ones(1,H);
end
log_c0 = log_c;

% max fun evals
%---------------
MaxFunEvals=options.MaxFunEvals;
if isempty(MaxFunEvals)
    MaxFunEvals=inf;
end

% max time
%---------------
MaxTime=options.MaxTime;
if isempty(MaxTime)
    MaxTime=inf;
end

% number of simulations
%-----------------------
N=options.N;
% NG=options.NG;

% target acceptance rate
%-------------------------
alpha=options.alpha;
alpha_swap = options.alpha_swap;

% thinning
%--------------
nthin = options.nthin;

% visualization/monitoring
%--------------------------
verbose_iter = options.verbose_iter;

%---------------------------------------------------------
rwm_fixed_p=options.rwm_fixed_p;
rwm_exp = options.rwm_exp;
sw_exp_rm = options.sw_exp_rm;
ram_adapt = options.ram_adapt;
separate_shape_adaptation = options.separate_shape_adaptation || ram_adapt;

% Initial temperature differences
%---------------------------------
Dlog_T = ones(1,H-1);% options.Dlog_T;% Initial log diff. of temperatures
lambda = inv_temperatures();

% LAMBDA=lambda(ones(N+1,1),:);
% The initial points for all levels
%-----------------------------------
if isempty(mu)
    mu=.5*(ub+lb);
end

if size(mu,2) == 1
    x = mu(:,ones(1,H));
else
    x = mu;
end

% Evaluate the initial values of the log-densities
fx = zeros(1,H);
funevals=0;
for istage=1:H
    fx(istage) = logf(x(:,istage));
    funevals=funevals+1;
    while is_violation(fx(istage),penalty)
        x(:,istage)=lb+(ub-lb).*rand(d,1);
        fx(istage) = logf(x(:,istage));
        funevals=funevals+1;
    end
    if ~isfinite(fx(istage))
        warning('Initial value of the log-density infinite')
    end
end

[bestx,bestf]=update_best(x,fx);
dd=@draw_parameter;
if isempty(SIG)
    nn=max(d+10,100);
    xsig=zeros(d,nn);
    fxsig=zeros(1,nn);
    fevals=zeros(1,nn);
    for icol=1:nn
        [xsig(:,icol),fxsig(icol),fevals(icol)]=dd();
    end
    funevals=funevals+sum(fevals);
    SIG=cov(xsig.');% SIG=1e-4*eye(d);
    [bestx,bestf,xx_,fxx_]=update_best([x,xsig],[fx,fxsig]);
    x=xx_(:,1:H);
    fx=fxx_(1:H);
end
Rchol=chol(SIG,'lower');

% Initial values for the mean, covariance and scalings
R0 = Rchol;
if separate_shape_adaptation
    m = x;
    if size(R0,3) == 1
        R0 = repmat(R0, [1, 1, H]);
    end
else
    m = mean(x,2);
end
Rchol = R0;

% Initialise some auxiliary variables
acceptance_prob = zeros(1,H);

% Array containing the (thinned) simulated samples
X = struct('f',{},'x',{});

% Acceptance statistics for RWM and switch:
accepted_metropolis_total = zeros(1,H);
accepted_swaps_total = zeros(1,H-1);
swaps_number_total = zeros(1,H-1);

% Auxiliary variables for the mean acceptance probability
% of switch between levels
N_sw = zeros(1,H-1);

metrop_func=@metropolis_draw;

% dating of events
%-----------------
t0=clock;

% amount of time spent
%---------------------
tic

isim=0;
iter=0;
not_converged=isim<N && funevals< MaxFunEvals && etime(clock,t0)<MaxTime;
while not_converged
    iter=iter+1;
    
    % gaussian proposal
    %---------------------
    u_randn = randn(d,H);
    
    % Metropolis moves
    %------------------
    fixed_component = (rand(1,H) <= rwm_fixed_p);
    stage_one = 1;
    fevals=zeros(1,H);
    for istage=1:H % parfor istage=1:H
        which_stage = separate_shape_adaptation*istage+...
            (1-separate_shape_adaptation)*stage_one;
        if fixed_component(istage)
            sqrt_cSIG = exp(log_c0(istage))*R0(:,:,which_stage)';
        else
            sqrt_cSIG = exp(log_c(istage))*Rchol(:,:,which_stage)';
        end
        [x(:,istage),fx(istage),...
            acceptance_prob(istage),accepted,u_randn(:,istage),fevals(istage)]=...
            metrop_func(x(:,istage),fx(istage),sqrt_cSIG,u_randn(:,istage),istage);
        accepted_metropolis_total(istage) = accepted_metropolis_total(istage) + accepted;
    end
    funevals=funevals+sum(fevals);
    
    % covariance and/or tune adaptation
    %----------------------------------
    alpha_diff = (~fixed_component).*(acceptance_prob-alpha);
    gamma_ = (iter+1)^(-rwm_exp);
    if ram_adapt
        Rchol = ram_adapt_shape(Rchol, gamma_, u_randn, alpha_diff);
    else
        % Adapt the scales of the level
        log_c = log_c + gamma_*alpha_diff;
        % Do the location-shape adaptation:
        [m, Rchol] = mean_cov_adapt(m, Rchol, x, gamma_);
    end
    
    % Swap states between temperatures
    %----------------------------------
    [x, fx, acc_sw_, alpha_sw_, who_swapped] = one_random_swap(x, fx);
    N_sw = N_sw + who_swapped;
    accepted_swaps_total = accepted_swaps_total + acc_sw_;
    
    % Adaptation of the temperature schedule
    %---------------------------------------
    lambda = inv_temperatures(alpha_sw_);
    
    % Record the statistics
    swaps_number_total = swaps_number_total + N_sw;
    N_sw(:) = 0;
    %     alpha_sw_est(:) = 0;
    
    if rem(iter,nthin) == 0
        isim=isim+1;
        X(isim) = struct('f',fx(1),'x',x(:,1));
        % Be sure to include only the stored guys... there is always a
        % chance that the best guy is lost in thinning but that is unlikely
        % here 
        if (minimization && fx(1)<bestf)||(~minimization && fx(1)>bestf)
            bestf=fx(1);
            bestx=x(:,1);
        end
    end

    if rem(isim,verbose_iter)==0
        t = toc;
        fprintf('iter: %0.0f, bestf: %0.4f, funevals: %0.0f, ETA %.0fs; accept. rate metrop %.2f%%, accept. rate swaps %.2f%%  \n', ...
            iter, bestf, funevals, t/isim*(N-isim), mean(accepted_metropolis_total./isim)*100, ...
            mean(accepted_swaps_total./swaps_number_total)*100 );
    end
    
    not_converged=isim<N && funevals< MaxFunEvals && etime(clock,t0)<MaxTime;
    
end
swaps_number_total = swaps_number_total + N_sw;

stats = struct('acceptance_rate_rwm', accepted_metropolis_total/N,...
    'acceptance_rate_swaps', accepted_swaps_total./swaps_number_total,...
    'time_elapsed', toc);

Results=struct('X',X,'bestf',bestf,'bestx',bestx,...
    'lambda',lambda,...
    'm',m,...
    'SIG',transpose(Rchol(:,:,1))*Rchol(:,:,1),...
    'c',exp(log_c),...
    'stats',stats,...
    'funevals',funevals);

    function [x,fx,accept_prob,accepted,u,fevals]=metropolis_draw(x0,f0,S,u,istage)
        % draw from the gaussian distribution centered at theta_star
        % with distribution mu_i, c_i*SIG_i. c_i is chosen so that the
        % acceptance rate is within a target range
        [x,fx,fevals,u]=draw_parameter(x0,S,u);
        [x,fx,accept_prob,accepted]=metropolis_selection(x0,f0,x,fx);
        % % %         % dynamic thinning
        % % %         %-----------------
        % % %         while rand>ps
        % % %             % do not store, redraw
        % % %             xx=x+sqrt_cSIG*randn(d,1);
        % % %             bad=xx<lb; xx(bad)=lb(bad);
        % % %             bad=xx>ub; xx(bad)=ub(bad);
        % % %             fxx=logf(xx); funevals=funevals+1;
        % % %             [x,fx,accept_prob,accepted]=metropolis_selection(x,fx,xx,fxx);
        % % %         end
        
        function [xx,fxx,accept_prob,success]=metropolis_selection(x0,f0,x,fx)
            % first select by the strength of violation, then by the level of
            % the function
            %         disp('deb selection not yet implemented!!!')
            
            accept_prob=acceptance_probability(lambda(istage),fx,f0);
            
            u_rand=rand;
            success=u_rand<accept_prob;
            if success
                xx=x;
                fxx=fx;
            else
                xx=x0;
                fxx=f0;
            end
        end
    end

    function [x,fx,fevals,u]=draw_parameter(x0,S,u)
            chol_style=nargin>0;
            x=draw_one();
            [fx,retcode]=logf(x);
            fevals=1;
            
            while retcode||is_violation(fx,penalty)
                u=[];
                x=draw_one();
                [fx,retcode]=logf(x);
                fevals=fevals+1;
                if debug && retcode
                    decipher(retcode)
                end
            end
            
            function x=draw_one()
                if chol_style
                    if isempty(u)
                        u=randn(d,1);
                    end
                    x=x0+S*u;
                    bad=x<lb; x(bad)=lb(bad);
                    bad=x>ub; x(bad)=ub(bad);
                else
                    x=lb+(ub-lb).*rand(d,1);
                end
            end
        end
        
    function lambda = inv_temperatures(alpha_sw_)
        if fixed_temperatures
            if geometric_tempering
                rho=sqrt(eps)^(1/(H-1));
                lambda=rho.^(0:H-1);
            else
                lambda=fliplr(linspace(sqrt(eps),1,H));
            end
        else
            %-------------------------------
            if nargin
                gamma_temper = ((H-1)/sum(who_swapped))*(iter+1)^(-sw_exp_rm);
                Dlog_T = Dlog_T + gamma_temper*who_swapped.*(alpha_sw_ - alpha_swap);
            end
            %-------------------------------
            T = ones(1,H);
            for k = 1:H-1
                T(k+1) = T(k) + exp(Dlog_T(k));
            end
            lambda = 1./T;
        end
    end

    function Rchol = ram_adapt_shape(Rchol, gamma, u, alpha_diff)
        % Robust AM adaptation for each of the levels
        
        % ram_adapt_shape.m
        %
        u_norm_sq = sum(u.^2,1);
        
        for istage_ = 1:H
            const = sqrt(min(0.9,d*gamma)*abs(alpha_diff(istage_))/u_norm_sq(istage_));
            if alpha_diff(istage_) < 0
                sign_='-';
            else
                sign_='+';
            end
            R_ = Rchol(:,:,istage_);
            Rchol(:,:,istage_)=cholupdate_lower(R_,const*R_*u(:,istage_),sign_);
        end
        
    end

    function [m, Rchol] = mean_cov_adapt(m, Rchol, x, gamma)
        % The mean & covariance adaptation
        
        % mean_cov_adapt.m
        %
        C_mean = zeros(d);
        m_mean = zeros(d,1);
        
        for istage_ = 1:H
            if separate_shape_adaptation
                dx = x(:,istage_) - m(:,istage_);
                % Mean update
                m(:,istage_) = m(:,istage_) + gamma*dx;
                % Covariance
                Rchol(:,:,istage_)=cholupdate_lower(sqrt(1-gamma)*Rchol(:,:,istage_),sqrt(gamma)*dx);
            else
                dx = x(:,istage_)-m;
                m_mean = m_mean + dx/H;
                C_mean = C_mean + dx*dx'/H;
            end
        end
        
        if ~separate_shape_adaptation
            % Mean
            m = (1-gamma)*m + gamma*m_mean;
            % Covariance
            Rchol = chol((1-gamma)*(Rchol*Rchol.') + gamma*C_mean,'lower');
        end
        
    end

    function Rlow=cholupdate_lower(Rlow,v,sign_)
        if nargin<3
            sign_='+';
        end
        Rlow = transpose(cholupdate(transpose(Rlow),v,sign_));
    end

    function [x, fx, acc_sw, alpha_swap, who_swapped] = one_random_swap(x, fx)
        
        % one_random_swap.m
        %
        who_swapped = zeros(1,H-1);
        alpha_swap = who_swapped;
        acc_sw = alpha_swap;
        
        % Pick random integer in [1,H-1]
        picked_stage = min(floor(rand*(H-1)+1),H-1);
        
        [acc_sw(picked_stage), alpha_swap(picked_stage)] = try_swap(picked_stage);
        
        who_swapped(picked_stage) = 1;
        
        function [acc_sw, acceptance_prob] = try_swap(chosen_stage)
            % Try to swap the levels istage <-> istage+1
            
            % try_swap.m
            %
            % The acceptance probability:
            delta_lambda=lambda(chosen_stage)-lambda(chosen_stage+1);
            acceptance_prob=acceptance_probability(delta_lambda,fx(chosen_stage+1),fx(chosen_stage));
            
            % An accept-reject for the swap
            if rand <= acceptance_prob
                % Swap states
                x_ = x(:,chosen_stage);
                fx_ = fx(chosen_stage);
                
                x(:,chosen_stage) = x(:,chosen_stage+1);
                fx(chosen_stage) = fx(chosen_stage+1);
                
                x(:,chosen_stage+1) = x_;
                fx(chosen_stage+1) = fx_;
                
                % Update statistics
                acc_sw = 1;
            else
                acc_sw = 0;
            end
            
        end
        
    end

    function prob=acceptance_probability(delta,lfx,lfx0)
        lfx0=delta*(lfx-lfx0);
        if minimization
            lfx0=-lfx0;
        end
        prob = min(1,exp(lfx0));
    end

    function [bestx,bestf,x,fx]=update_best(x,fx)
        if minimization
            [fx,tags] = sort(fx,2);
        else
            [fx,tags] = sort(fx,2,'descend');
        end
        x=x(:,tags);
        bestx=x(:,1);
        bestf=fx(1);
    end
end


function flag=is_violation(f,penalty)
flag=abs(f)>=penalty||~isreal(f);
if ~isreal(f)
    warning('f not real!!!')
end
end

