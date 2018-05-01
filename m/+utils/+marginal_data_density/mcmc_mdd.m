function [log_mdd,extras] = mcmc_mdd(theta_draws,lb,ub,options)
% MCMC_MDD -- computes various types of log marginal data density
%
% ::
%
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws)
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws,lb)
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws,lb,ub)
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws,lb,ub,options)
%
% Args:
%
%    - **theta_draws** [struct]: with fields "f" and "x". Each parameter is
%    defined as a structure, which means that theta_draws is a vector of
%    structures. "x" is the parameter vector and "f" is the NEGATIVE of the
%    log posterior kernel evaluated at "x". In case "f" is instead the log
%    posterior kernel itself, option **maximization** below has to be set to
%    "true".
%
%    - **lb** [empty|vector]: lower bound of the search space. Necessary only
%    for the swz algorithm. Conveniently replaced with the lower bounds of
%    theta_draws if empty.
%
%    - **ub** [empty|vector]: upper bound of the search space. Necessary only
%    for the swz algorithm. Conveniently replaced with the upper bounds of
%    theta_draws if empty.
%
%    - **options** [struct]: with possible fields
%      - **log_post_kern** [function handle]: function computing the log
%      posterior kernel for a given parameter vector
%      - **center_at_mean** [{false}|true]: if true, the distribution is
%      centered at the mean. Else, it is centered at the mode, which should be
%      the maximum of the log posterior kernel in theta_draws
%      - **algorithm** [{mhm}|swz|mueller|bridge|is|ris|cj|laplace|laplace_mcmc]:
%          - **mhm** is the modified harmonic mean
%          - **swz** is the Sims, Waggoner and Zha (2008) algorithm
%          - **mueller** is the unpublished Mueller algorithm (see Liu,
%          Waggoner and Zha 2011).
%          - **bridge** is the Meng and Wong (1996) algorithm.
%          - **is** is the Importance sampling algorithm.
%          - **ris** is the reciprocal importance sampling algorithm.
%          - **cj** is the Chib and Jeliazkov (2001) algorithm.
%          - **laplace** is the laplace approximation
%          - **laplace_mcmc** is the laplace approximation using the
%          covariance of the mcmc draws rather than the one obtain from
%          computing and inverting the numerical hessian.
%      - **L** [{500}|integer|struct]: number of IID draws or IID draws
%      - **maximization** [{false}|true]: Informs the procedure about whether
%      we have a maximization or a minimization problem.
%      - **debug** [{false}|true]: print useful information during estimation.
%      - **mhm_tau** [{(.1:.1:.9)}|vector|scalar]: truncation probabilities
%      for the MHM algorithm
%      - **swz_pvalue** [{90}|scalar]: scalar for the computation of the lower
%      bound in the SWZ algorithm
%      - **bridge_TolFun** [numeric|{sqrt(eps)}]: convergence criterion in the
%      BRIDGE algorithm
%
% Returns:
%    :
%
%    - **log_mdd** [numeric]: log marginal data density
%
%    - **extras** [empty|struct]: further output from specific algorithms,
%    which will include the iid draws as they can be re-use in some other
%    algorithm.
%
% Note:
%
% Example:
%
%    See also:

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
algorithms={'mhm','swz','mueller','bridge','is','ris','cj','laplace','laplace_mcmc'};
defaults={ % arg_names -- defaults -- checks -- error_msg
    'log_post_kern',[],@(x)isa(x,'function_handle'),'log_posterior_kern should be a function handle'
    'center_at_mean',false,@(x)isscalar(x) && islogical(x),'center_at_mean should be a logical scalar'
    'algorithm','mhm',@(x)any(strcmp(x,algorithms)),['algorithm must be one of ',cell2mat(strcat(algorithms,'|'))]
    'L',500,@(x)num_fin_int(x)||isstruct(x),'L (# of IID draws) should be an integer or a structure containing draws'
    'maximization',false,@(x)isscalar(x) && islogical(x),'maximization should be a logical scalar'
    'debug',false,@(x)isscalar(x) && islogical(x),'debug should be a logical scalar'
    'mhm_tau',(.1:.1:.9),@(x)isempty(x)||(all(x>0) && all(x<1)),'mhm_tau (truncation probabilities of MHM) must all be in (0,1)'
    'swz_pvalue',90,@(x)num_fin(x) && x>=0 && x<=100,'swz_pvalue (???) must be in [0,100]'
    'bridge_TolFun',sqrt(eps),@(x)num_fin(x) && x>0,'bridge_TolFun (tolerance level) should be a positive scalar'
    'bridge_fix_point',true,@(x)isscalar(x) && islogical(x),'bridge_fix_point should be a logical scalar'
    'bridge_initialize_is',false,@(x)isscalar(x) && islogical(x),'bridge_initialize_is should be a logical scalar'
    };

if nargin<4
    options=[];
    if nargin<3
        ub=[];
        if nargin<2
            lb=[];
            if nargin<1
                log_mdd=cell2struct(defaults(:,2),defaults(:,1),1);
                return
            end
        end
    end
end

if isempty(options)
    options=struct();
end

options=parse_arguments(defaults,options);

center_at_mean=options.center_at_mean;

log_post_kern=options.log_post_kern;

algorithm=options.algorithm;

debug=options.debug;

bridge_TolFun=options.bridge_TolFun;

bridge_fix_point=options.bridge_fix_point;

bridge_initialize_is=options.bridge_initialize_is;

is_need_log_posterior_kernel=~strcmp(algorithm,'mhm');

if isempty(log_post_kern) && is_need_log_posterior_kernel
    error(['log_post_kern cannot be empty for algorithm "',algorithm,'"'])
end

% By default we do minimization : check whether the user wants something
% else
%------------------------------------------------------------------------
maximization=options.maximization;
thecoef=2*maximization-1;
% multiply accordingly
LogPost_M=thecoef*[theta_draws.f];
[~,best_loc]=max(LogPost_M);

% the exponential explodes very quickly. We apply the following fix
%-------------------------------------------------------------------
facility=true;
facilitator=facility*LogPost_M(best_loc);

theta_mode=theta_draws(best_loc(1));

theta_draws=[theta_draws.x];

[d,M] = size(theta_draws);

if isempty(lb)
    lb=min(theta_draws,[],2);
end
if ~isequal(size(lb),[d,1])
    error('wrong format of lower bound')
end

if isempty(ub)
    ub=max(theta_draws,[],2);
end
if ~isequal(size(ub),[d,1])
    error('wrong format of upper bound')
end

% recentering function
%----------------------
recenter=@(x)utils.optim.recenter(x,lb,ub,3);

% get moments
%--------------
[x_bar,SigmaMode,Sigma_i,det_Sigma]=moments();
Shat = chol(SigmaMode,'lower');

NumWorkers=utils.parallel.get_number_of_workers();
extras=[];
L = options.L;
iid_draws_=[];
if isstruct(L)
    iid_draws_=L;
    L=numel(L);
end
switch algorithm
    case 'swz'
        log_mdd=do_sims_waggoner_zha_2008();
    case 'mhm'
        log_mdd=do_modified_harmonic_mean();
    case 'mueller'
        log_mdd=do_mueller_2004();
    case 'bridge'
        log_mdd=do_bridge_1996();
    case 'is'
        log_mdd=do_importance_sampling();
    case 'ris'
        log_mdd=do_reciprocal_importance_sampling();
    case 'cj'
        log_mdd=do_chib_jeliazkov_2001();
    case 'laplace'
        log_mdd=do_laplace();
    case 'laplace_mcmc'
        log_mdd=do_laplace_mcmc();
end

extras.iid_draws=iid_draws_;

% Methods
%----------

    function log_mdd=do_laplace()
        H=utils.hessian.finite_differences(log_post_kern,theta_mode.x);
        if maximization
            H=-H;
        end
        Hinv=H\eye(size(H,1));
        log_mdd=utils.marginal_data_density.laplace_mdd(...
            LogPost_M(best_loc),Hinv);
    end

    function log_mdd=do_laplace_mcmc()
        Hinv=cov(theta_draws.');
        log_mdd=utils.marginal_data_density.laplace_mdd(...
            LogPost_M(best_loc),Hinv);
    end

    function log_mdd=do_chib_jeliazkov_2001()
        % Posterior kernel at the mode
        %--------------------------------------------------
        numerator=LogPost_M(best_loc);
        
        % weighting-function density for the M original draws
        %-----------------------------------------------------
        Logq_M=old_draws_in_weighting_function('CJ');
        
        % log posterior kernel of the IID draws
        %---------------------------------------
        [LogPost_L]=iid_draws('CJ');
        
        % acceptance probability for the IID draws relative to the
        % posterior mode 
        %-----------------------------------------------------------------
        alpha=zeros(1,L);
        for idraw=1:L
            logr=LogPost_L(idraw)-numerator;
            alpha(idraw)=exp(min(0,logr));
        end
        
        extras=struct('cj_alpha',alpha);
        
        % approximate log posterior at the mode
        %--------------------------------------
        big=max(Logq_M);
        top=exp(Logq_M-big);
        lp_thet_y=big+log(mean(top))-log(mean(alpha));
        
        % finally compute the marginal data density
        %------------------------------------------
        log_mdd=numerator-lp_thet_y;
    end

    function log_mdd=do_reciprocal_importance_sampling()
        % Sylvia Frühwirth-Schnatter (2004): "Estimating marginal
        % likelihoods for mixture and Markov switching models using bridge
        % sampling techniques". Econometrics Journal 7(1) pp 143--167
        Logq_M=old_draws_in_weighting_function('RIS');
        ratio = Logq_M-LogPost_M;
        ratiomax=max(ratio);
        log_mdd = -1*(ratiomax+log(mean(exp(ratio-ratiomax))));
    end

    function log_mdd=do_importance_sampling(LogPost_L,Logq_L)
        % Sylvia Frühwirth-Schnatter (2004): "Estimating marginal
        % likelihoods for mixture and Markov switching models using bridge
        % sampling techniques". Econometrics Journal 7(1) pp 143--167
        if nargin<2
            [LogPost_L,Logq_L]=iid_draws('IS');
        end
        ratio = LogPost_L-Logq_L;
        ratiomax=max(ratio);
        log_mdd = ratiomax+log(mean(exp(ratio-ratiomax)));
    end

    function log_mdd=do_bridge_1996()
        % Xiao-Li Meng and Wing Hung Wong (1996): " Simulating Ratios of
        % Normalizing Constants via a Simple Identity: A Theoretical
        % Exploration". Statistica Sinica, 6, 831860.
        
        % 1- draw L vectors and evaluate them at the posterior kernel and
        % at the weighting function
        %------------------------------------------------------------------
        [LogPost_L,Logq_L]=iid_draws('BRIDGE');
        
        % 2- for the M original draws compute the density from the
        % weighting function
        %----------------------------------------------------------
        Logq_M=old_draws_in_weighting_function('BRIDGE');
        
        lmdd0=0;
        if bridge_initialize_is
            lmdd0=do_importance_sampling(LogPost_L,Logq_L);
        end
        
        if bridge_fix_point
            [log_mdd,lmdd,conv]=fix_point_strategy(lmdd0);
        else
            [log_mdd,lmdd,conv]=iterative_strategy(lmdd0);
        end
        extras=struct('bridge_lmdd',lmdd,'bridge_conv',conv);

        function [log_mdd,lmdd,conv]=iterative_strategy(lmdd0)
            % 3- for a number of iterations, update an initial guess for the
            % log MDD
            %----------------------------------------------------------------
            max_iter=1000;
            lmdd=zeros(max_iter,1);
            lmdd(1)=lmdd0;
            conv=inf(max_iter,1);
            iter=1;
            while conv(iter)>bridge_TolFun
                iter=iter+1;
                if iter==max_iter
                    lmdd(iter+1000)=0;
                    conv(iter+1000)=inf;
                    max_iter=max_iter+1000;
                end
                lmdd(iter)=do_one_iteration(lmdd(iter-1));
                conv(iter)=abs(lmdd(iter)-lmdd(iter-1));
                if debug
                    fprintf(1,'BRIDGE: iter #: %0.0f, lmdd: %0.4f, lmdd: %0.4f\n',...
                        iter,lmdd(iter),conv(iter));
                end
                if iter==2
                    % use a reasonable number at the start
                    conv(iter-1)=conv(iter);
                end
            end
            log_mdd=lmdd(iter);
            lmdd=lmdd(1:iter);
            conv=conv(1:iter);
        end
        
        function [log_mdd,lmdd,conv]=fix_point_strategy(lmdd0)
            fzero_options=optimset('Display','none');
            if debug
                fprintf('BRIDGE: Now bisecting...');
                tic
                fzero_options=optimset('Display','iter');
            end
            next=0;
            max_iter=1000;
            conv=zeros(max_iter,1);
            lmdd=zeros(max_iter,1);
            
            log_mdd=fzero(@one_right_hand_side,lmdd0,fzero_options);
            if debug
                fprintf('Done in %0.4d seconds\n',toc);
            end
            conv=conv(1:next);
            lmdd=lmdd(1:next);
            
            function rhs=one_right_hand_side(lmdd_iter)
                % numerator
                %-----------
                rq=alpha_times_p_or_q(Logq_L,LogPost_L-lmdd_iter,2);
                
                % denominator
                %-------------
                rmc=alpha_times_p_or_q(Logq_M,LogPost_M-lmdd_iter,1);
                
                % numerator minus denominator
                %----------------------------
                max_rq=max(rq);
                max_rmc=max(rmc);
                
                rhs=max_rq+log(mean(exp(rq-max_rq)))-...
                    (max_rmc+log(mean(exp(rmc-max_rmc))));
                next=next+1;
                if next==max_iter
                    conv(next+1000)=0;
                    lmdd(next+1000)=0;
                    max_iter=max_iter+1000;
                end
                conv(next)=rhs;
                lmdd(next)=lmdd_iter;
            end
        end
        
        function a_p=alpha_times_p_or_q(f1,f2,stud)
            biggest=max([f1;f2],[],1);
            f1=f1-biggest;
            f2=f2-biggest;
            if stud==1
                a_p=f1;
            else
                a_p=f2;
            end
            a_p=a_p-log(...
                exp(log(L)+f1)+...
                exp(log(M)+f2)...
                );
        end
        
        function lmdd=do_one_iteration(lmdd)
            % numerator
            %-----------
            lp_q=LogPost_L-lmdd;
            rq=alpha_times_p_or_q(Logq_L,lp_q,2);
            
            % denominator
            %-------------
            lp_m=LogPost_M-lmdd;
            rmc=alpha_times_p_or_q(Logq_M,lp_m,1);
            
            % numerator minus denominator
            %----------------------------
            max_rq=max(rq);
            max_rmc=max(rmc);
            
            lmdd=lmdd+...
                max_rq+log(mean(exp(max(rq-max_rq,-1e15))))-...
                (max_rmc+log(mean(exp(max(rmc-max_rmc,-1e15)))));
        end
    end

    function log_mdd=do_sims_waggoner_zha_2008()
        swz_pvalue=options.swz_pvalue;
        % step 1: compute the r(s)
        % ------------------------
        r=compute_r();
        
        % step 2: compute the parameters of the h-density
        % ------------------------------------------------
        [aa,bb,v]=compute_a_b_v();
        
        % step 3: draw L parameters from the elliptical distribution
        % -----------------------------------------------------------
        % the following quantity will be re-used many times
        bv_av=bb^v-aa^v;
        
        [~,~,Lswz,qLswz]=elliptical_distribution_log_posterior_kernel_and_bounds();
        extras=struct('swz_L',Lswz,'swz_qL',qLswz);
        
        % step 4: Log values of the elliptical distribution at the posterior draws
        % -------------------------------------------------------------------------
        log_Integration_constant = gammaln(d/2)-(log(2)+d/2*log(pi)+...
            1/2*log(det(SigmaMode)));
        
        logG_theta = log_Integration_constant+log(v/bv_av)+(v-d)*log(r);
        
        % step 5: Compute the log marginal data density
        % ----------------------------------------------
        logh = -log(qLswz)+logG_theta(:).';
        
        % finally compute the modified harmonic mean
        %--------------------------------------------
        good=LogPost_M>Lswz & isfinite(logh);
        log_mdd = sum(exp(logh(good)-LogPost_M(good)+facilitator))/M;
        
        log_mdd = facilitator-log(log_mdd);
        
        function [lpost_ell,theta_ell,Lswz,qLswz]=elliptical_distribution_log_posterior_kernel_and_bounds()
            lpost_ell = zeros(1,L);
            v_i=1/v;
            norm2=@(x)sqrt(sum(x.^2));
            theta_ell=zeros(d,L);
            % we parallelize the following because of the expensive
            % function call log_post_kern
            if debug
                fprintf('SWZ: Now sampling and evaluating elliptical draws...');
                tic
            end
            parfor (ii=1:L,NumWorkers)
                u=rand;
                rd=(bv_av*u)^v_i;
                x = randn(d,1);
                draw=rd/norm2(x)*Shat*x+x_bar;
                theta_ell(:,ii)=recenter(draw);
                % compute the posterior kernel values for the iid draws
                lpost_ell(ii) = thecoef*log_post_kern(theta_ell(:,ii)); %#ok<PFBNS>
            end
            if debug
                fprintf('Done in %0.4d seconds\n',toc);
            end
            % lower bound from the posterior kernel based such that
            % prob(f<=Lswz)=swz_pvalue
            % -------------------------------------------------------------
            [Lswz,qLswz]=find_lower_bound_and_quantile();
            function [Lswz,qLswz]=find_lower_bound_and_quantile()
                SortLogPost = sort(LogPost_M);
                Lswz = SortLogPost(ceil((100-swz_pvalue)/100*M));%
                qLswz = sum(lpost_ell>Lswz)/L;
            end
        end
        
        function [a,b,v]=compute_a_b_v()
            rDist = sort(r);
            c1 = rDist(ceil((1/100)*M));
            c10 = rDist(ceil((1/10)*M));
            c90 = rDist(floor((9/10)*M));
            %
            v = log(1/9)/log(c10/c90);
            b = c90/(0.9^(1/v));
            a = c1;
        end
        
        function r=compute_r()
            r= zeros(M,1);
            for ii=1:M;
                dev = theta_draws(:,ii)-x_bar;
                r(ii) = sqrt(dev'*Sigma_i*dev);
            end
        end
    end

    function log_mdd=do_mueller_2004()
        % Logq_L: weighting draws evaluated through weighting function
        %--------------------------------------------------------------
        % LogPost_L: weighting draws evaluated through the log-posterior kernel
        %----------------------------------------------------------------
        [LogPost_L,Logq_L]=iid_draws('MUELLER');
        
        % old draws evaluated through the weighting function
        %----------------------------------------------------
        Logq_M=old_draws_in_weighting_function('MUELLER');

        % get the ratios
        %----------------
        log_rh=LogPost_L-Logq_L;
        log_rp=Logq_M-LogPost_M;
        
        fzero_options=optimset('Display','none');
        if debug
            fprintf('MUELLER: Now bisecting...');
            tic
            fzero_options=optimset('Display','iter');
        end
        next=0;
        max_iter=1000;
        record=zeros(max_iter,1);
        log_mdd=fzero(@hfunc,0,fzero_options);
        if debug
            fprintf('Done in %0.4d seconds\n',toc);
        end
        extras=struct('mueller_convergence',record(1:next));
        
        log_mdd=thecoef*log_mdd;
        
        function h=hfunc(log_c)
            % left side
            %------------
            left=log_c+log_rh;
            good=left<0;
            left=sum(1-exp(left(good)))/L;
            % right side
            %------------
            right=log_rp-log_c;
            good=right<0;
            right=sum(1-exp(right(good)))/M;
            h=left-right;
            next=next+1;
            if next==max_iter
                record(next+1000)=0;
                max_iter=max_iter+1000;
            end
            record(next)=h;
        end
    end

    function log_mhm=do_modified_harmonic_mean()
        mhm_tau=options.mhm_tau;
        critical_value=chi2inv(mhm_tau,d);
        log_tau=log(mhm_tau);
        post_tmn=zeros(1,numel(log_tau));
        % we serialize the following in order to save on memory and because
        % calling conditional_likelihood should be relatively cheap
        if debug
            fprintf('MHM: Just doing it...');
            tic
        end
        for idraw=1:M % parfor (idraw=1:M,NumWorkers)
            % evaluate proposal
            %------------------
            [loglik,v_iF_v]=normal_weighting(theta_draws(:,idraw)-x_bar);
            % geweke's truncated multivariate normal
            %---------------------------------------
            good=v_iF_v<=critical_value;
            log_f=loglik-log_tau(good);
            post_tmn(good)=post_tmn(good)+exp(log_f-LogPost_M(idraw)+facilitator);
        end
        if debug
            fprintf('Done in %0.4d seconds\n',toc);
        end
        post_tmn=post_tmn/M;
        log_mhm=max(facilitator-log(post_tmn));
    end

% Utility functions
%--------------------
    function [x_bar,Sigma_bar,Sigma_i,det_Sigma]=moments()
        if debug
            fprintf('MDD: Now computing moments from %0.0d draws...',M);
            tic
        end
        if center_at_mean
            x_bar = (1/M)*sum(theta_draws,2);
        else
            x_bar = theta_mode.x;
        end
        Sigma_bar = 1/M*(theta_draws*theta_draws')-(x_bar*x_bar');
        Sigma_bar=utils.cov.project(Sigma_bar);
        Sigma_i=Sigma_bar\eye(d);
        det_Sigma=det(Sigma_bar);
        if debug
            fprintf('Done in %0.4d seconds\n\n',toc);
        end
    end

    function [loglik,v_iF_v,lik]=normal_weighting(v)
        v_iF_v=v'*Sigma_i*v;
        loglik=-.5*(d*log(2*pi)+log(det_Sigma)+v_iF_v);
        if nargout>2
            lik=exp(loglik);
        end
    end

    function [LogPost_L,Logq_L]=iid_draws(caller)
        if debug
            fprintf('%s: Now evaluating %0.0d IID draws...',caller,L);
            tic
        end
        eval_weighting=nargout>1;
        if eval_weighting
            Logq_L=zeros(1,L);
        end
        is_drawn=~isempty(iid_draws_);
        for idraw=1:L
            if ~is_drawn
                if idraw==1
                    iid_draws_=struct();
                end
                v=Shat*randn(d,1);
                draw=recenter(x_bar+v);
                iid_draws_(idraw).x=draw;
                iid_draws_(idraw).f=log_post_kern(draw);
            end
            if eval_weighting
                v=iid_draws_(idraw).x-x_bar;
                Logq_L(idraw)=normal_weighting(v);
            end
        end
        LogPost_L=thecoef*[iid_draws_.f];
        if debug
            fprintf('Done in %0.4d seconds\n\n',toc);
        end
    end

    function Logq_M=old_draws_in_weighting_function(caller)
        if debug
            fprintf('%s: Now evaluating %0.0d old draws in weighting function...',caller,M);
            tic
        end
        Logq_M=zeros(1,M);
        for idraw=1:M
            draw=theta_draws(:,idraw);
            v=draw-x_bar;
            Logq_M(idraw)=normal_weighting(v);
        end
        if debug
            fprintf('Done in %0.4d seconds\n\n',toc);
        end
    end
end

