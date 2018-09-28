function [log_mdd,extras]=swz(obj,swz_pvalue,opts)

n=nargin;

set_inputs()

moments(obj,opts)

m=obj.moms.m;

theta_draws=obj.theta_draws;

M=obj.M;

d=obj.d;

LogPost_M=obj.LogPost_M;

facilitator=obj.facilitator;

% step 1: compute the r(s)
% ------------------------
r=compute_r();

% step 2: compute the parameters of the h-density
% ------------------------------------------------
[aa,bb,v]=compute_a_b_v();

% step 3: draw opts.L parameters from the elliptical distribution
% -----------------------------------------------------------
% the following quantity will be re-used many times
bv_av=bb^v-aa^v;

[~,~,Lswz,qLswz]=elliptical_distribution_log_posterior_kernel_and_bounds();

extras=struct('swz_L',Lswz,'swz_qL',qLswz);

% step 4: Log values of the elliptical distribution at the posterior draws
% -------------------------------------------------------------------------
log_Integration_constant = gammaln(d/2)-(log(2)+d/2*log(pi)+...
    1/2*log(det(m.Sigma_bar)));

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
        
        lpost_ell = zeros(1,opts.L);
        
        v_i=1/v;
        
        norm2=@(x)sqrt(sum(x.^2));
        
        theta_ell=zeros(d,opts.L);
        
        % we parallelize the following because of the expensive
        % function call log_post_kern
        if opts.debug
            
            fprintf('SWZ: Now sampling and evaluating elliptical draws...\n');
            
            tic
            
        end
        
        NumWorkers=utils.parallel.get_number_of_workers();
        
        parfor (ii=1:opts.L,NumWorkers)
            
            u=rand;
            
            rd=(bv_av*u)^v_i;
            
            x = randn(d,1);
            
            draw=rd/norm2(x)*m.Shat*x+m.x_bar;
            
            theta_ell(:,ii)=obj.recenter(draw); %#ok<*PFBNS>
            
            % compute the posterior kernel values for the iid draws
            lpost_ell(ii) = obj.thecoef*obj.log_post_kern(theta_ell(:,ii)); 
            
        end
        
        if opts.debug
            
            fprintf('Done in %0.4f seconds\n',toc);
            
        end
        
        % lower bound from the posterior kernel based such that
        % prob(f<=Lswz)=swz_pvalue
        % -------------------------------------------------------------
        [Lswz,qLswz]=find_lower_bound_and_quantile();
        
        function [Lswz,qLswz]=find_lower_bound_and_quantile()
            
            SortLogPost = sort(obj.LogPost_M);
            
            Lswz = SortLogPost(ceil((100-swz_pvalue)/100*M));%
            
            qLswz = sum(lpost_ell>Lswz)/opts.L;
            
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
        
        for ii=1:M
            
            dev = theta_draws(:,ii)-m.x_bar;
            
            r(ii) = sqrt(dev'*m.Sigma_i*dev);
            
        end
        
    end

    function set_inputs()
        
        if n<3
            
            opts=[];
            
            if n<2
                
                swz_pvalue=[];
                
            end
            
        end
        
        if isempty(swz_pvalue)
            
            swz_pvalue=90;
            
        end
        
        num_fin=@(x)isnumeric(x) && isscalar(x);
        
        test=@(x)num_fin(x) && x>=0 && x<=100;
        
        assert(test(swz_pvalue),'swz_pvalue (???) must be in [0,100]')
        
        opts=mdd.global_options(opts);
        
    end

end
