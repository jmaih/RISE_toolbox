function [log_mdd,qLswz,Lswz] = mcmc_mdd(theta_draws,lb,ub,options)
% MCMC_MDD -- computes various types of log marginal data density
%
% Syntax
% -------
% ::
%
%   [log_mdd,qLswz,Lswz] = MCMC_MDD(theta_draws)
%
%   [log_mdd,qLswz,Lswz] = MCMC_MDD(theta_draws,lb)
%
%   [log_mdd,qLswz,Lswz] = MCMC_MDD(theta_draws,lb,ub)
%
%   [log_mdd,qLswz,Lswz] = MCMC_MDD(theta_draws,lb,ub,options)
%
% Inputs
% -------
%
% - **theta_draws** [struct]: with fields "f" and "x". Each parameter is
% defined as a structure, which means that theta_draws is a vector of
% structures. "x" is the parameter vector and "f" is the NEGATIVE of the
% log posterior kernel evaluated at "x". In case "f" is instead the log
% posterior kernel itself, option **maximization** below has to be set to
% "true".
%
% - **lb** [empty|vector]: lower bound of the search space. Necessary only
% for the swz algorithm. Conveniently replaced with the lower bounds of
% theta_draws if empty.
%
% - **ub** [empty|vector]: upper bound of the search space. Necessary only
% for the swz algorithm. Conveniently replaced with the upper bounds of
% theta_draws if empty.
%
% - **options** [struct]: with possible fields
%   - **log_post_kern** [function handle]: function computing the log
%   posterior kernel for a given parameter vector
%   - **center_at_mean** [{false}|true]: if true, the distribution is
%   centered at the mean. Else, it is centered at the mode, which should be
%   the maximum of the log posterior kernel in theta_draws
%   - **algorithm** [{mhm}|swz|mueller|bridge]: 
%       - **mhm** is the modified harmonic mean
%       - **swz** is the Sims, Waggoner and Zha (2008) algorithm
%       - **mueller** is the unpublished Mueller algorithm (see Liu,
%       Waggoner and Zha 2011). 
%       - **bridge** is the Meng and Wong (1996) algorithm. 
%   - **mhm_tau** [{(.1:.1:.9)}|vector|scalar]: truncation probabilities for
%   the modified harmonic mean 
%   - **L** [{[]}|integer]: number of elliptical draws (active if swz is
%   true) 
%   - **swz_pvalue** [{90}|scalar]: scalar for the computation of the lower
%   bound in the SWZ algorithm
%   - **maximization** [{false}|true]: Informs the procedure about whether
%   we have a maximization or a minimization problem.
%   - **debug** [{false}|true]: print useful information during estimation.
%
% Outputs
% --------
%
% - **log_mdd** [numeric]: log marginal data density
%
% - **qLswz** [numeric]: with an estimate of the fraction of iid draws from
% the elliptical distribution such that the posterior kernel values at
% these draws are greater than Lswz. (empty if swz is false)
%
% - **Lswz** [numeric]: such that swz_pvalue percent of the posterior draws have a
% log posterior kernel values greater then Lswz. (empty if swz is false)
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
algorithms={'mhm','swz','mueller','bridge','is','ris'};
defaults={ % arg_names -- defaults -- checks -- error_msg
    'swz_pvalue',90,@(x)num_fin(x) && x>=0 && x<=100,'swz_pvalue (???) must be in [0,100]'
    'center_at_mean',false,@(x)isscalar(x) && islogical(x),'center_at_mean should be a logical scalar'
    'log_post_kern',[],@(x)isa(x,'function_handle'),'log_posterior_kern should be a function handle'
    'L',500,@(x)num_fin_int(x),'L (# of IID draws) should be an integer'
    'mhm_tau',(.1:.1:.9),@(x)isempty(x)||(all(x>0) && all(x<1)),'mhm_tau (truncation probabilities of MHM) must all be in (0,1)'
    'maximization',false,@(x)isscalar(x) && islogical(x),'maximization should be a logical scalar'
    'algorithm','mhm',@(x)any(strcmp(x,algorithms)),['algorithm must be one of ',cell2mat(strcat(algorithms,'|'))]
    'debug',false,@(x)isscalar(x) && islogical(x),'debug should be a logical scalar'
    'bridge_maxiter',1000,@(x)num_fin_int(x),'bridge_maxiter number of bridge iterations'
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

options=utils.miscellaneous.parse_arguments(defaults,options);

center_at_mean=options.center_at_mean;

log_post_kern=options.log_post_kern;

algorithm=options.algorithm;

debug=options.debug;

bridge_maxiter=options.bridge_maxiter;

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
Lswz=[];qLswz=[];
L = options.L;
switch algorithm
    case 'swz'
        [log_mdd,Lswz,qLswz]=do_sims_waggoner_zha_2008();
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

    function log_mdd=do_importance_sampling()
        % Sylvia Frühwirth-Schnatter (2004): "Estimating marginal
        % likelihoods for mixture and Markov switching models using bridge
        % sampling techniques". Econometrics Journal 7(1) pp 143--167
        [LogPost_L,Logq_L]=iid_draws('IS');
        ratio = LogPost_L-Logq_L;
        ratiomax=max(ratio);
        log_mdd = ratiomax+log(mean(exp(ratio-ratiomax)));
    end

    function log_mdd=do_bridge_1996()
        % Xiao-Li Meng and Wing Hung Wong (1996): " Simulating Ratios of
        % Normalizing Constants via a Simple Identity: A Theoretical
        % Exploration". Statistica Sinica, 6, 831–860.
        
        % 1- draw L vectors and evaluate them at the posterior kernel and
        % at the weighting function
        %------------------------------------------------------------------
        [LogPost_L,Logq_L]=iid_draws('BRIDGE');
        
        % 2- for the M original draws compute the density from the
        % weighting function
        %----------------------------------------------------------
        Logq_M=old_draws_in_weighting_function('BRIDGE');
        
        % 3- for a number of iterations, update an initial guess for the
        % log MDD
        %----------------------------------------------------------------
        lmdd=zeros(bridge_maxiter,1);
        for iter=2:bridge_maxiter
            lmdd(iter)=do_one_iteration(lmdd(iter-1));
            if debug
                fprintf(1,'BRIDGE: iter #: %0.0f, lmdd: %0.4f\n',...
                    iter,lmdd(iter));
            end
        end
        
        log_mdd=lmdd(end);
        
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
        end
    end

    function [log_mdd,Lswz,qLswz]=do_sims_waggoner_zha_2008()
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
        log_mdd=fzero(@hfunc,0,fzero_options);
        if debug
            fprintf('Done in %0.4d seconds\n',toc);
        end
        
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
        LogPost_L=zeros(1,L);
        Logq_L=zeros(1,L);
        for idraw=1:L
            v=Shat*randn(d,1);
            draw=recenter(x_bar+v);
            LogPost_L(idraw)=thecoef*log_post_kern(draw);
            Logq_L(idraw)=normal_weighting(v);
        end
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

