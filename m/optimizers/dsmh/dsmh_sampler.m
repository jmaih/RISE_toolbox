function [Params,ff,log_mdd]=dsmh_sampler(logf,lb,ub,options,mu,SIG)
% try and save the draws at each stage in a way that simulation can proceed
% with exactly the same settings if interrupted... basically in a way that
% one can just load an continue. This means the tempering levels, etc will
% be stored too.
%
% The optimization in this function is a maximization and not a
% minimization as typically done... of course that can be put as an option.

% Reference:
% Daniel F. Waggoner, Hongwei Wu and Tao Zha (2014): "Dynamic Striated
% Metropolis-Hastings Sampler for High-Dimensional Models"
%
% - We modify the initial distribution to be the multivariate uniform
% distribution rather than the multivariate t-student as suggested by the
% authors.

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
typical_nshocks=20;
typical_sample=200;
defaults={ % arg_names -- defaults -- checks -- error_msg
    'nu',5,@(x)num_fin_int(x) && x>0,' nu (degrees of freedom of multivariate t-student) must be a positive integer'
    'M',20,@(x)num_fin_int(x) && x>=1,'M (# striations) must be a positive integer'
    'minimization',true,@(x)isscalar(x) && islogical(x),'minimization must be a logical scalar'
    'c',1,@(x)num_fin(x) && x>0,'c (tuning parameter) should be a positive scalar'
%     'recombine_all',false,@(x)isscalar(x) && islogical(x),'recombine_all should be a logical scalar'
    'NG',5000,@(x)num_fin(x) && x>0,'NG(# simulations) should be a strictly positive integer'
    'G',1,@(x)num_fin(x) && x>0,'G(# cores) should be a strictly positive integer'
    'H',40,@(x)num_fin(x) && x>0,'H(# tempering stages) should be a strictly positive integer'
    'p',[],@(x)isempty(x)||(num_fin(x) && x>=0 && x<1),'p (prob of drawing from prev. distrib.) must be in [0,1) when not empty. If empty, it will be set to 0.1*ps'
    'lambda_1',1/(10*typical_nshocks*typical_sample),@(x)num_fin(x) && x>0 && x<=1,'lambda_1 (first level of tempering) should be in (0,1]'
    'geometric_lambda',true,@(x)isscalar(x) && islogical(x),'geometric_lambda should be a logical scalar'
    'ps',.9,@(x)num_fin(x) && x>0 && x<=1,'ps (prob of storing a metropolis draw) must be in (0,1]'
    'alpha',.234,@(x)num_fin(x) && x>0 && x<=1,'alpha (target acceptance rate) must be in (0,1]'
    };
if nargin==0
    Params=cell2struct(defaults(:,2),defaults(:,1),1);
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

% option for whether the optimization is a maximization or a minimization
% ------------------------------------------------------------------------
minimization=options.minimization;

% degrees of freedom of multivariate t-student distribution
%-----------------------------------------------------------
nu=options.nu;

% number of striations
%----------------------
M=options.M;

% ess_iw_min=0.1;

% tuning parameter for the covariance of the metropolis
%-------------------------------------------------------
c=options.c;

% % recombination of particiles/bees/vectors
% %-------------------------------------------
% % if recombine_all==false, we leave each striation unchanged... the only
% % communication with the other striations is through the common covariance
% % matrix
% recombine_all=options.recombine_all;

% number of simulations
%-----------------------
NG=options.NG;

% number of cores
%-----------------
% default should be 1
G=options.G;

% number of stages
%------------------
H=options.H;

% probability of storing a metropolis draw
%-------------------------------------------
ps=options.ps;

% probability of drawing from the distribution of the previous stage.
%----------------------------------------------------------------------
% Instead of just drawing from the previous stage, one can devise a
% bee-type of algorithm within the striation and still use the metropolis
% to accept or reject
p=options.p;
if isempty(p)
    p=0.1*ps;
end

% target acceptance rate
%-------------------------
alpha=options.alpha;

% Tempering levels
%------------------
geometric_lambda=options.geometric_lambda;

lambda_vector=tempering_steps(options.lambda_1);

% draws per striation
%---------------------
% choose L so that the there is an equal number of draws in each striation
ndraws_per_striation=ceil(NG/M);

NG=ndraws_per_striation*M;

ndraws_per_striation_per_core=floor(ndraws_per_striation/G);
remains=ndraws_per_striation-ndraws_per_striation_per_core*G;
ndraws_per_striation_per_core=ndraws_per_striation_per_core*ones(1,G);
ndraws_per_striation_per_core(1:remains)=ndraws_per_striation_per_core(1:remains)+1;

% draw initial distribution: striations in rows, cores in columns
%-----------------------------------------------------------------
[STRIATIONS_CORES_X,STRIATIONS_CORES_F,log_I_old]=initial_distribution(true);

sqrt_cSIG=[];

for istage=1:H
    
    lambda_old=lambda_vector(istage);
    
    lambda=lambda_vector(istage+1);
    
    [w,wtilde]=importance_weights();
    
    % compute log marginal data density
    %-----------------------------------
    log_I=log_I_old+log(mean(wtilde));
    
    % update moments to be used in the metropolis step.
    %----------------------------------------------------
    update_moments(w);
    
    % draw theta, ensuring that an equal number of draws lies in each
    % striation
    NEW_X=STRIATIONS_CORES_X;
    NEW_F=STRIATIONS_CORES_F;
    METROP_NDRAWS=cell(M,G);
    METROP_NACCEPT=cell(M,G);
    for icore=1:G
        ndraws_i=ndraws_per_striation_per_core(icore);
        if minimization
            bestf=inf;
        else
            bestf=-inf;
        end
        for st=1:M
            n_metrop_draws=0;
            n_accept=0;
            [theta_star,f_star]=first_mean_metropolis(st,icore);
            for idraw_=1:ndraws_i
                if rand<=p
                    % draw randomly from previous distribution: the drawn
                    % theta has to belong to the same striation as theta_star
                    [NEW_X{st,icore}(:,idraw_),NEW_F{st,icore}(idraw_)]=...
                        striation_draw(st,icore,idraw_);
                else
                    % make a metropolis draw
                    [NEW_X{st,icore}(:,idraw_),NEW_F{st,icore}(idraw_),success]=...
                        metropolis_draw(theta_star,f_star);
                    theta_star=NEW_X{st,icore}(:,idraw_);
                    f_star=NEW_F{st,icore}(idraw_);
                    n_metrop_draws=n_metrop_draws+1;
                    n_accept=n_accept+success;
                end
            end
            if minimization
                bestf=min(bestf,NEW_F{st,icore}(idraw_));
            else
                bestf=max(bestf,NEW_F{st,icore}(idraw_));
            end
            METROP_NDRAWS{st,icore}=n_metrop_draws;
            METROP_NACCEPT{st,icore}=n_accept;
        end
        % adjust the scale within the core? If so do it here
        % after a certain number of iterations
        fprintf('core %0.0f, stage %0.0f of %0.0f, lambda %8.4f, best f(x) %8.4f \n\n',...
            icore,istage,H,lambda,bestf);
    end
    % adjust the scale using all cores? If so do it here. WWZ seem to have
    % a unique ci for each generation
    
    % gather results
    %----------------
    fthet=cell2mat(NEW_F);
    n_metrop_draws=sum(sum(cell2mat(METROP_NDRAWS)));
    n_accept=sum(sum(cell2mat(METROP_NACCEPT)));
    
    if minimization
        bestf=min(fthet(:));
    else
        bestf=max(fthet(:));
    end
    
    fprintf('stage %0.0f of %0.0f, lambda %8.4f, best f(x) %8.4f \n\n',istage,H,lambda,bestf);
    
    % next round
    %------------
    log_I_old=log_I;
    STRIATIONS_CORES_X=NEW_X;
    STRIATIONS_CORES_F=NEW_F;
end
% collect the results
%---------------------
offset=0;
Params=zeros(d,NG);
ff=zeros(1,NG);
for st=1:M
    for icore=1:G
        np=numel(STRIATIONS_CORES_F{st,icore});
        ff(offset+(1:np))=STRIATIONS_CORES_F{st,icore};
        Params(:,offset+(1:np))=STRIATIONS_CORES_X{st,icore};
        offset=offset+np;
    end
end

log_mdd=log_I;

    function [theta_star,f_star]=first_mean_metropolis(st,icore)
        theta_star=mean(STRIATIONS_CORES_X{st,icore},2);
        f_star=logf(theta_star);
%         [theta_star,f_star]=striation_draw(st,icore);
    end

    function [x,fx]=striation_draw(st,icore,ii)
%         pos=(st-1)*ndraws_per_striation+randi(ndraws_per_striation);
        x=generate_mutant(STRIATIONS_CORES_X{st,icore},ii);
        fx=logf(x);
    end

    function mutant=generate_mutant(xx,ii)
        nx=size(xx,2);
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

    function [x,fx,accepted]=metropolis_draw(x0,f0)
        % draw from the gaussian distribution centered at theta_star
        % with distribution mu_i, c_i*SIG_i. c_i is chosen so that the
        % acceptance rate is within a target range
        x=x0+sqrt_cSIG*randn(d,1);
        bad=x<lb; x(bad)=lb(bad);
        bad=x>ub; x(bad)=ub(bad);
        fx=logf(x);
        [x,fx,accepted]=metropolis_selection(x0,f0,x,fx);
        % dynamic thinning
        %-----------------
        while rand>ps
            % do not store, redraw
            xx=x+sqrt_cSIG*randn(d,1);
            bad=xx<lb; xx(bad)=lb(bad);
            bad=xx>ub; xx(bad)=ub(bad);
            fxx=logf(xx);
            [x,fx,accepted]=metropolis_selection(x,fx,xx,fxx);
        end
    end

    function [xx,fxx,success]=metropolis_selection(x0,f0,x,fx)
        % first select by the strength of violation, then by the level of
        % the function
%         disp('deb selection not yet implemented!!!')
        lfx0=lambda*(fx-f0);
        if minimization
            lfx0=-lfx0;
        end
        accept_prob=min(1,exp(lfx0));
        success=rand<accept_prob;
        if success
            xx=x;
            fxx=fx;
        else
            xx=x0;
            fxx=f0;
        end
    end

%     function ess=effective_sample_size(w)%ESS_IW
%         ess=1/sum(w.^2);
%     end

    function [w,wtilde]=importance_weights()
        fthet_=cell2mat(STRIATIONS_CORES_F);
        wtilde=(lambda-lambda_old)*fthet_(:);
        wtilde=exp(wtilde);
        w=wtilde/sum(wtilde);
        % core-specific quantities
        %--------------------------
    end

    function update_moments(w)
        mu=0;
        SIG=0;
        ii=0;
        for icore_=1:G
            for iv=1:ndraws_per_striation_per_core(icore_)
                % for each column, run through all the striations
                %-------------------------------------------------
                for istri=1:M
                    theta_i=STRIATIONS_CORES_X{istri,icore_}(:,iv);
                    ii=ii+1;%(icore_-1)*M+iv
                    wthet=w(ii)*theta_i;
                    mu=mu+wthet;
                    SIG=SIG+wthet*theta_i.';
                end
            end
        end
        SIG=SIG-(mu*mu.');
        % we may need to procrustify this...
        SIG=utils.cov.nearest(SIG);
        
        update_scaled_choleski();
        
        function update_scaled_choleski
            if istage>1
                acceptance_rate=n_accept/n_metrop_draws;
                cnew=utils.mcmc.retune_coefficient(c,[alpha,alpha],acceptance_rate);
                if cnew==0
                    keyboard
                end
            end
            sqrt_cSIG=chol(c*SIG,'lower');
        end
    end

    function [CORE_PARTIONS_X,CORE_PARTIONS_F,log_I_init]=initial_distribution(tstudent)
        if nargin==0
            tstudent=false;
        end
        % this has to be a probability density and not a Kernel and so it
        % sums to 1.
        I_init=1;
        log_I_init=log(I_init);
        if tstudent
            if isempty(SIG)
                SIG=eye(d);
            end
            if isempty(mu)
                mu=.5*(ub+lb);
            end
            iSIG=SIG\eye(d);
            logDetSIG=log(det(SIG));
            A=chol(SIG,'lower');
            % draws from the multivariate t-student distribution. See
            % Gelman, Carlin, Stern and Rubin (2004, page 581)
            sqrt_nu_over_x=sqrt(nu./chi2rnd(nu,1,NG));
            theta__=mu(:,ones(1,NG))+A*randn(d,NG).*sqrt_nu_over_x(ones(d,1),:);
            fthet__=zeros(1,NG);
            for ivec=1:NG
                fthet__(ivec)=multivariate_tstudent_density(theta__(:,ivec));
            end
        else
            ub_minus_lb=ub-lb;
            theta__=bsxfun(@plus,lb,bsxfun(@times,ub_minus_lb,rand(d,NG)));
            % log density of the multivariate uniform distribution
            fthet__=-log(prod(ub_minus_lb));
            if isnan(fthet__)||isinf(fthet__)
                % do the t-student option, which may require some
                % iterations that might be implemented later or not...
                [CORE_PARTIONS_X,CORE_PARTIONS_F,log_I_init]=initial_distribution(true);
                return
            else
                fthet__=fthet__(1,ones(1,NG));
            end
            % no need to sort since they all have the same height... Note
            % those draws are not necessarily feasible with respect to the
            % objective function, which is a potential problem.
        end
        % same number of guys in each striation
        [~,tags]=sort(fthet__);
        badges=ones(ndraws_per_striation,1)*(1:M);
        badges=badges(:).';
        newtags(tags)=1:numel(tags);
        badges=badges(newtags);
        CORE_PARTIONS_X=cell(M,G);
        CORE_PARTIONS_F=cell(M,G);
        for it=1:M
            partition_i=badges==it;
            badges(partition_i)=[];
            
            CORE_PARTIONS_X(it,:)=mat2cell(theta__(:,partition_i),d,ndraws_per_striation_per_core); 
            theta__(:,partition_i)=[];
            
            CORE_PARTIONS_F(it,:)=mat2cell(fthet__(partition_i),1,ndraws_per_striation_per_core); 
            fthet__(partition_i)=[];
        end
                
        function [lnp]=multivariate_tstudent_density(theta__)
            thet_mu=theta__-mu;
            lnp=gammaln(.5*(nu+1))-(gammaln(.5*nu)+.5*d*log(nu+pi))...
                -0.5*logDetSIG-0.5*(nu+d)*log(1+1/nu*thet_mu.'*iSIG*thet_mu);
        end
    end

    function lambda_vector=tempering_steps(lambda_1)
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
    end
end
