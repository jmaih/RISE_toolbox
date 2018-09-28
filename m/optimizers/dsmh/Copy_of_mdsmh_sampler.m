function [Results]=mdsmh_sampler(logf,lb,ub,options,mu,SIG)
% Matlab complains that there is no documentation. But then here is one!!!
% To do:
% - beggar-thy-neighbor: whenever you sample a new guy, you put him in the
% striation it belongs i.e. in the striation with the same level of energy
% - striation specific metropolis, i.e. each striation has its own
% metropolis point, initialized at each iteration to be the mean of the
% guys in the striation: this breaks correlation
% - re-introduce geometric weights such that the guys with low energy have
% a wider striation?
% - we could also consider randomizing the weights or having weights depend
% on the performance?
% - finish the equi-energy sampler, doing it with structures like in this
% function. check the adaptive equi-energy sampler for tips for
% improvement.
% - finish the parallel tempering sampler and optimize it, again write in
% as structures. check for the latest developments
% - re-write the metropolis hastings algorithm with the same philosophy as
% in this function and make it independent of rise objects
% - re-write the marginal data densities: 
%   - Meng and Wong
%   - Adaptive Importance Sampling,
%   - Mueller, 
%   - Chib-Jeliazkov, 
%   - Waggoner and Zha, 
%   - Geweke
% - Include moves from the metaheuristic literature in the sampler? I most
% likely would be blamed for insulting convergence theorems and
% reversibility of markov processes?
% - Do a paper that
%   - presents the problem
%   - scale intervals if the product goes to infinity?
%   - argues that the techniques are also and perhaps more useful for
%   nonlinear models estimated using frequentist approaches
%   - runs examples on two models
%       - the Schorfheide model
%       - the Liu-Waggoner-Zha model. Try to see whether the simulator
%       recover the posterior mode
%   - plot both the bivariate marginal and univariate marginal and argue
%   that they could be multimodal, while the metropolis would want to have
%   them unimodal.
% - Remove field posterior simulation ? where to put the posterior
% statistics?
% - Does it have to be the case that the number of draws should hold in
% memory for all the exercises including the calculation of MDD, IRFS,
% intervals, etc. or does it have to be like dynare where we collect the
% info from the disk as necessary?
% - simulators independent
% - put metropolis in the simulators department alongside EE, DSMH, etc.
% - parfor with condition all intensive processes:
%   - first detect the presence of slaves
%   - curvature
%   - derivatives?
%   - plots
%   - scale 
%   - curvature

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
typical_nshocks=20;
typical_sample=200;
defaults={ % arg_names -- defaults -- checks -- error_msg
%     'nu',5,@(x)num_fin_int(x) && x>0,' nu (degrees of freedom of multivariate t-student) must be a positive integer'
    'M',20,@(x)num_fin_int(x) && x>=1,'M (# striations) must be a positive integer'
    'minimization',true,@(x)isscalar(x) && islogical(x),'minimization must be a logical scalar'
    'c',1,@(x)num_fin(x) && x>0,'c (tuning parameter) should be a positive scalar'
    'fixed_scaling',false,@(x)isscalar(x) && islogical(x),'fixed_scaling (fixed tuning parameter) should be a logical scalar'
    'N',5000,@(x)num_fin(x) && x>0,'N(# simulations per core) should be a strictly positive integer'
    'G',1,@(x)num_fin(x) && x>0,'G(# cores) should be a strictly positive integer'
    'H',40,@(x)isempty(x)||(num_fin(x) && x>0),'H(# tempering stages) should be a strictly positive integer or empty'
    'p',[],@(x)isempty(x)||(num_fin(x) && x>=0 && x<1),'p (prob of drawing from prev. distrib.) must be in [0,1) when not empty. If empty, it will be set to 0.1*ps'
    'lambda_1',1/(10*typical_nshocks*typical_sample),@(x)num_fin(x) && x>0 && x<=1,'lambda_1 (first level of tempering) should be in (0,1]'
    'geometric_lambda',true,@(x)isscalar(x) && islogical(x),'geometric_lambda should be a logical scalar'
    'ps',.9,@(x)num_fin(x) && x>0 && x<=1,'ps (prob of storing a metropolis draw) must be in (0,1]'
    'alpha',.234,@(x)num_fin(x) && x>0 && x<=1,'alpha (target acceptance rate) must be in (0,1]'
%     'initial_distribution','uniform',@(x)any(strcmpi(x,{'t-student','uniform','laplace'})),'initial_distribution must be "t-student","uniform" or "laplace"'
    'common_moments',false,@(x)isscalar(x) && islogical(x),'common_moments (common moments across cores) should be a logical scalar'
    'rwm_exp',0.6,@(x)num_fin(x) && x>=0,'rwm_exp (Exponent of random-walk adaptation step size) must be a positive scalar'
    'ess_min',0.1,@(x)num_fin(x) && x>0 && x<1, 'ess_min (constant for controlling the evolution of lambda) should be in (0,1)'
    'common_lambda',true,@(x)isscalar(x) && islogical(x),'common_lambdas (common lambda across cores) should be a logical scalar'
    'p_mutant',0,@(x)num_fin(x) && x>=0 && x<=1,'p_mutant (prob of drawing a mutant) must be in [0,1]'
    'energy_adjust',2,@(x)num_fin_int(x) && x>0,'energy_adjust (adjustment factor for striation bounds) must be a positive number'
    'penalty',[],@(x)isempty(x)||num_fin(x) && x>0,'penalty (worst value possible in absolute value) must be empty or a finite positive number'
    };

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

if isempty(mu)
    mu=.5*(ub+lb);
end

if isempty(SIG)
    SIG=eye(d);
end

options=utils.miscellaneous.parse_arguments(defaults,options);

if ischar(logf)
    logf=str2func(logf);
end

if ~isa(logf,'function_handle')
    error('logf should be a function handle or a string')
end

% worst value that the function can assume
%------------------------------------------
if isempty(options.penalty)
    if options.minimization
        options.penalty=1e+6;
    else
        options.penalty=-1e+6;
    end
end
penalty=options.penalty;

% number of striations
%----------------------
M=options.M;

% tuning parameter for the covariance of the metropolis
%-------------------------------------------------------
log_c=log(options.c);

% number of cores
%-----------------
% default should be 1
G=options.G;

% number of simulations
%-----------------------
N=options.N;

% draws per striation
%---------------------
% choose L so that the there is an equal number of draws in each striation
ndraws_per_striation=ceil(N/M);

N=ndraws_per_striation*M;

% draw initial distribution: 
%----------------------------

[genpop,striations,striation_ranges,log_I_old]=initial_draws();

% Tempering levels
%------------------
geometric_lambda=options.geometric_lambda;

% number of stages
%------------------
H=options.H;
fixed_horizon=~isempty(H);
if fixed_horizon
    options.common_lambda=true;
end

% common lambda across cores?
%------------------------------
common_lambda=options.common_lambda;

% common moments across cores?
%------------------------------
common_moments=options.common_moments;

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

% non-standard options to run on the cores: do this after the options have
% been updated
%--------------------------------------------------------------------------
options_dispatch=dispatch_options();

log_I=cell(1,G+1);

istage=0;
done=false;
while ~done
    istage=istage+1;
    
    lambda_old=lambda_vector(istage);
    
    lambda=lambda_vector(istage+1);
    % do each core
    %--------------------
    parfor (icore=1:G,G*(G>1)) % use parallel here only if G>1
        if common_moments
            % use the common mu and... in that case we may want to say that
            % they are not to be computed inside each core...
            options_dispatch{icore}.SIG=SIG;
            options_dispatch{icore}.mu=mu;
        end
        if fixed_horizon||common_lambda
            options_dispatch{icore}.lambda(istage+1)=lambda;
        elseif lambda_converged(options_dispatch{icore}.lambda(end))
            % If we were done in the previous stage, no need to proceed
            % and the lambda is computed within each core
            continue
        end
        options_dispatch{icore}.istage = istage;
        [genpop{icore},log_I{icore+1},options_dispatch{icore}]=do_one_core(...
            logf,genpop{icore},log_I_old{icore+1},options_dispatch{icore}); %#ok<PFOUS>
    end
    
    F=cell2mat(genpop); F=[F.f];
    
    [w,wtilde,ess]=importance_weights(F,lambda,lambda_old,options.minimization);
    
    % compute log marginal data density across all cores
    %----------------------------------------------------
    log_I{1}=log_I_old{1}+log(mean(wtilde(:)));
    
    % update moments to be used in the metropolis step.
    %----------------------------------------------------
    [mu,SIG]=update_moments(genpop,w);
    
    log_I_old=log_I;
    
    if options.minimization
        bestf=min(matricize(F));
    else
        bestf=max(matricize(F));
    end
    
    metrop_draws=0;
    metrop_accept=0;
    funevals=0;
    for icore=1:G
        metrop_draws=metrop_draws+options_dispatch{icore}.metrop_draws;
        metrop_accept=metrop_accept+options_dispatch{icore}.metrop_accept;
        funevals=funevals+options_dispatch{icore}.funevals;
    end
    accept_ratio=metrop_accept/metrop_draws;
    display_progress([],istage,options.H,lambda,funevals,bestf,...
        accept_ratio,ess,numel(w))
    
    if fixed_horizon
        done=istage==H;
        % the lambda is already computed
    else
        if common_lambda
            % the horizon is not fixed but the lambda decision is still common:
            % all cores will converge at the same time.
            done=lambda_converged(lambda);
        else
            % each core is responsible for computing its own lambda: the
            % different cores most likely will converge at different horizons.
            % We are done only when all cores have converged
            %------------------------------------------------------------------
            done=true;
            for icore=1:G
                done=done && lambda_converged(options_dispatch{icore}.lambda(end));
                if ~done
                    break
                end
            end
        end
        if ~done
            % update the grand lambda at least for the computation of overall
            % moments even if the individual cores use specific moments
            %-----------------------------------------------------------------
            lambda_vector(istage+2)=find_next_lambda(lambda,matricize(F),...
                options.minimization,options.ess_min);
        end
    end
end

% remove the common options in each suboption
%---------------------------------------------
fields=fieldnames(options);
for icore=1:G
    if icore==1
        trimmed_options=rmfield(options_dispatch{icore},fields);
    else
        trimmed_options(icore)=rmfield(options_dispatch{icore},fields);
    end
end
genpop=matricize(genpop);
F=[genpop.f];
icol=find(F==bestf,1,'first');
Results=struct('X',genpop,...
    'FX',F,'bestf',bestf,...
    'bestx',genpop(:,icol).x,...
    'lambda',lambda_vector,...
    'm',mu,...
    'SIG',SIG,...
    'funevals',funevals,...
    'log_mdd',matricize(log_I),...
    'ess',ess,...
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
        else
            lambda_vector=[0,lambda_1];
        end
    end

    function [bigpop,striations,striation_ranges,log_I_init]=initial_draws()
        
        % this has to be a probability density and not a Kernel and so it
        % sums to 1. Do initialization for each core and for the general...
        I_init=1;
        log_I_init=repmat({log(I_init)},1,G+1);
        
        % Total number of simulations
        NG=N*G;
        
        ub_minus_lb=ub-lb;
        theta=bsxfun(@plus,lb,bsxfun(@times,ub_minus_lb,rand(d,NG)));
        genpop_=cell(1,NG);
        parfor idraw=1:NG
            fx=logf(theta(:,idraw)); %#ok<PFBNS>
            while is_violation(fx,penalty)
                theta(:,idraw)=lb+ub_minus_lb.*rand(d,1);
                fx=logf(theta(:,idraw));
            end
            genpop_{idraw}=struct('f',fx,'x',theta(:,idraw));
        end
        ldens=-sum(log(ub_minus_lb));%ldens=log(1/prod(ub_minus_lb));
        if options.minimization
            ldens=-ldens;
        end
        genpop_=cell2mat(genpop_);
            
        % same number of guys in each striation
        [~,tags]=sort_population(genpop_);
        badges=ones(ndraws_per_striation*G,1)*(1:M);
        badges=badges(:).';
        newtags(tags)=1:numel(tags);
        badges=badges(newtags);
        bigpop=cell(1,G);
        striations=create_striations([genpop_.f],M,options.energy_adjust);
        striation_ranges=cell(1,M);
        for it=1:M
            striation_ranges{it}=(it-1)*ndraws_per_striation+1:it*ndraws_per_striation;
            partition_i=badges==it;
            badges(partition_i)=[];
            pop_part=genpop_(partition_i);
            genpop_(partition_i)=[];
            offset=0;
            for ig=1:G
                this_range=offset+(1:ndraws_per_striation);
                bigpop{ig}=[bigpop{ig},pop_part(this_range)];
                offset=offset+ndraws_per_striation;
                if it==M
                    % change the f to satisfy the requirements of a density
                    % function
                    [bigpop{ig}.f]=deal(ldens);
                end
            end
        end
    end

    function opt=dispatch_options()
        % add non-standard options to run on the cores. First the common
        %----------------------------------------------------------------
        opt=options;
        opt.SIG=SIG;
        opt.log_c=log_c;
        opt.lb=lb;
        opt.ub=ub;
        opt.funevals=0;
        opt.lambda=lambda_vector;
        opt.striations=striations(:).';
        opt.striation_ranges=striation_ranges;
        % now for the specifics
        %----------------------
        opt=repmat({opt},1,G);
        for ig=1:G
            theta_star=genpop{ig}(1:0);
            for istri=1:M
                lucky=randi(ndraws_per_striation);
                pos=striation_ranges{istri}(lucky);
                theta_star(1,istri)=genpop{ig}(pos);
            end
            opt{ig}.theta_star=theta_star;
            opt{ig}.core=ig;
        end
    end

end

function [pop,tags]=sort_population(pop)
[~,tags]=sort([pop.f]);
pop=pop(tags);
end

function flag=is_violation(f,penalty)
flag=abs(f)>=penalty||~isreal(f);
if ~isreal(f)
    warning('f not real!!!')
end
end

function [st,striations,isnew]=find_striation(fx,striations,F,energy_adjust)
if fx<striations(end) && fx>=striations(1)
    st=find(striations<=fx,1,'last');
    isnew=false;
else
    isnew=true;
    M=numel(striations)-1;
    if fx<striations(1)
        F(1)=fx;
    else
        F(end)=fx;
    end
    striations=create_striations(F,M,energy_adjust);
    st=find_striation(fx,striations);
end
end

function striations=create_striations(F,M,energy_adjust)

% each striation is such that x in [lb,ub) that is, lb <= x < ub in other
% words, if x ==ub, then it does not belong to the striation [lb,ub)

F=sort(F(:).');

nf=numel(F);

b=linspace(0,nf,M+1)+1;

if any(b-fix(b)~=0)
    error('unbalanced number of F and M')
end

striations=[F(b(1:end-1)),F(end)+energy_adjust];

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

function display_progress(core,istage,H,lambda,funevals,bestf,accept_ratio,ess,N,nsecs)
if nargin<10
    nsecs=nan;
end
string=['stage %0.0f of %0.0f, lambda %8.4f, funevals %4.0f, best f(x) %8.4f'...
    ', acceptance rate %4.2f, ESS(N) %4.4f(%0.0f), this iteration %0.4f sec'];
if isempty(H)
    H=nan;
end
data={istage,H,lambda,funevals,bestf,100*accept_ratio,ess,N,nsecs};
if isempty(core)
    string=['Aggregate:: ',string];
else
    string=['core %0.0f, ',string];
    data=[{core},data];
end
fprintf([string,'\n\n'],data{:});

end

function Y=matricize(Y)
if iscell(Y)
    Y=cell2mat(Y);
end
end

function [lamb,fval,exitflag,output,jacobian]=find_next_lambda(lambda_old,F,minimization,ess_min,solver)
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
            [~,~,ess]=importance_weights(F,lambda,lambda_old,minimization);
            discrep=ess-N_ess_min;
        end
    end
end

function [w,wtilde,ess]=importance_weights(F,lambda,lambda_old,minimization)
% the weights may collapse if the lambda differential is so wide that in
% case F has large components, the exponential goes to infinity. This
% probably happens because the initial distribution is not a probability
% density...
%-----------------------------------------------------------------------
F=matricize(F);
if minimization
    F=-F;
end
wtilde=(lambda-lambda_old)*F;
% log_big=max(wtilde);
wtilde=exp(wtilde);% wtilde=exp(wtilde-log_big);
w=wtilde/sum(wtilde);
if nargout>2
    % effective sample size
    %------------------------
    ess=effective_sample_size(w);
end

end

function ess=effective_sample_size(w)
ess=1/sum(w(:).^2);
end

function [mu,SIG]=update_moments(pop,w)

pop=matricize(pop);

mu=0;
SIG=0;
[~,ndraws]=size(pop);
for iv=1:ndraws
    theta_i=pop(iv).x;
    wthet=w(iv)*theta_i;
    mu=mu+wthet;
    SIG=SIG+wthet*theta_i.';
end
SIG=SIG-(mu*mu.');

% we may need to procrustify this...
SIG=utils.cov.nearest(SIG);
end

function log_c=update_scaling(log_c,fixed_scaling,gam,accept_ratio,alpha_range)

alpha_1=min(alpha_range);

alpha_2=max(alpha_range);

alpha=.5*(alpha_1+alpha_2);

alpha_diff = (~fixed_scaling).*(accept_ratio-alpha);

log_c = log_c + gam*alpha_diff;

end

function [pop,log_I,options]=do_one_core(logf,pop0,log_I_old,options)
p=options.p;
ps=options.ps;
probability_of_not_storing_metropolis=1-ps;
p_mutant=options.p_mutant;
minimization=options.minimization;
theta_star=options.theta_star(1);
lb=options.lb;
ub=options.ub;
istage=options.istage;
lambda_old=options.lambda(istage);
lambda=options.lambda(istage+1);
striations=options.striations;
striation_ranges=options.striation_ranges;
M=options.M; % =numel(striations)-1;
penalty=options.penalty;

in_funevals=0;

% separate in striations
%-------------------------
M_times_ndraws_per_striation=numel(pop0);
d=size(lb,1);
ndraws_per_striation=M_times_ndraws_per_striation/M;
pop=pop0;


sqrt_cSIG=chol(exp(options.log_c)*options.SIG,'lower');

% sqrt_cSIG=chol(exp(options.log_c)*options.SIG,'lower');
metrop_draws=0;
metrop_accept=0;
tic
for idraw=1:M_times_ndraws_per_striation
    theta=dsmh_draw();
    
    % replace the guy in the same position
    %--------------------------------------
    
    pop(idraw)=theta;
        
%     % find the right striation for the newcomer
%     [target_stri,striations,isnew]=find_striation(fx,striations,F,options.energy_adjust);
%     
%     if isnew
%         % striations have changed and so some vectors may now belong to
%         % a different striation. In that case we need to sort the
%         % database
%         [F,tags]=sort(F);
%         X=X(:,tags);
%     end
    
%     % place the draw in the striation it belongs. But then since there
%     % is room for only a limited number of draws per striation, we need
%     % to pop one guy and this could be:
%     % - the guy with the highest correlation (with the guy entering the
%     % database) among those with lower energy...
%     % - the guy with the lowest energy... this is the option we
%     % implement for now
%     
%     Ftarg=F(striation_ranges{target_stri});
%     if minimization
%         nailed=find(Ftarg==max(Ftarg),1,'first');
%     else
%         nailed=find(Ftarg==min(Ftarg),1,'first');
%     end
%     nailed=striation_ranges{target_stri}(nailed);
%     F(nailed)=fx;
%     X(:,nailed)=x;
end

% resort the database
%---------------------
pop=sort_population(pop);

% update for next round
%-----------------------
options.theta_star=theta_star;

% refresh the striations
%-------------------------
options.striations=striations;

[w,wtilde,options.ess]=importance_weights([pop.f],lambda,lambda_old,minimization);

log_I=log_I_old+log(mean(wtilde(:)));

accept_ratio=nan;
if metrop_draws
    gam_log_c = (istage+1)^(-options.rwm_exp);
    accept_ratio=metrop_accept/metrop_draws;
    options.log_c=update_scaling(options.log_c,options.fixed_scaling,...
        gam_log_c,accept_ratio,options.alpha);
end

options.metrop_draws=metrop_draws;
options.metrop_accept=metrop_accept;

if minimization
    bestf=pop(1).f;
else
    bestf=pop(end).f;
end

time_it_took=toc;
options.funevals=options.funevals+in_funevals;
display_progress(options.core,istage,options.H,lambda,...
    options.funevals,bestf,accept_ratio,options.ess,numel(w),time_it_took)

% update moments to be used in the metropolis step. Update these moments
% only if they are not common across cores. The mean will not be used
% directly, but the covariance will. The advantage of using a common
% covariance is that the dispersion is potentially more important.
%-------------------------------------------------------------------------
if ~options.common_moments
    [options.mu,options.SIG]=update_moments(pop,w);
end
% % % % if accept_ratio==0||any(isnan(options.SIG(:)))||any(options.SIG(:)>1e+4)
% % % %     save bad_things_happened
% % % %     error('no accepted draws')
% % % % end

% find the next lambda iff the lambda is not common across cores
%----------------------------------------------------------------
if ~options.common_lambda && ~lambda_converged(lambda)
    options.lambda(istage+2)=find_next_lambda(lambda,[pop.f],minimization,options.ess_min);
end

    function [theta]=dsmh_draw()
        
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
            if is_violation(theta.f,penalty)
                theta=theta0;
                accepted=false;
            else
                lfratio=lambda*(theta.f-theta0.f);
                lfg_ratio=lfratio+lgratio;
                if minimization
                    lfg_ratio=-lfg_ratio;
                end
                % Metropolis criterion
                %-----------------------
                crit=exp(min(0,lfg_ratio));
                accepted=rand<crit;
                if accepted
                    if is_metrop_draw
                        metrop_accept=metrop_accept+1;
                    end
                else
                    theta=theta0;
                end
            end
        end
        
        function [theta,lgratio]=previous_level_draw()
            nail=randi(ndraws_per_striation);
            % find the right striation for the newcomer
            [target_stri,striations]=find_striation(theta_star.f,...
                striations,[pop.f],options.energy_adjust);
            
            if rand<p_mutant
                [theta]=draw_mutant();
            else
                this_guy=striation_ranges{target_stri}(nail);
                theta=pop0(this_guy);
            end
            lgratio=lambda_old*(theta_star.f-theta.f);
            
            function [theta]=draw_mutant()
                xd=generate_mutant([pop0(striation_ranges{target_stri}).x],nail,lb,ub);
                theta.f=logf(xd);
                theta.x=xd;
                in_funevals=in_funevals+1;
            end
        end
        
        function [theta,lgratio]=gaussian_draw()
            xd=theta_star.x+sqrt_cSIG*randn(d,1);
            bad=xd<lb; xd(bad)=lb(bad);
            bad=xd>ub; xd(bad)=ub(bad);
            theta.f=logf(xd);
            theta.x=xd;
            in_funevals=in_funevals+1;
            lgratio=0;
            metrop_draws=metrop_draws+1;
        end
    end
end