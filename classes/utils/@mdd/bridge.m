function [log_mdd,extras]=bridge(obj,fix_point,opts)

% Xiao-Li Meng and Wing Hung Wong (1996): " Simulating Ratios of
% Normalizing Constants via a Simple Identity: A Theoretical
% Exploration". Statistica Sinica, 6, 831860.

n=nargin;

set_inputs()

iid_draws_=[];

nold=obj.M;

c=obj.thecoef;

if opts.draws_are_iid
    
    iid_draws_=struct();
    
    for idraw=1:opts.L
        
        iid_draws_(idraw).x=obj.theta_draws(:,end-idraw+1);
        
        % If maximization, c=1. If minimization, c=-1 and the current
        % operation will be reversed in iid_draws;
        iid_draws_(idraw).f=c*obj.LogPost_M(end-idraw+1);
        
    end
    
    nold=nold-opts.L;
    
end
% 1- draw opts.L vectors and evaluate them at the posterior kernel and
% at the weighting function
%------------------------------------------------------------------
[LogPost_L,Logq_L]=iid_draws(obj,iid_draws_,opts,'BRIDGE');

% 2- for the M original draws compute the density from the
% weighting function
%----------------------------------------------------------
Logq_M=old_draws_in_weighting_function(obj,opts,'BRIDGE');

Logq_M=Logq_M(1:nold);

LogPost_M=obj.LogPost_M(1:nold);

% initialize with importance sampling...
%---------------------------------------
lmdd0=0;% lmdd0=is(obj,[],opts);

max_iter=opts.fix_point_maxiter;

bridge_TolFun=opts.fix_point_TolFun;

lmdd=zeros(max_iter,1);

conv=inf(max_iter,1);

if fix_point
    
    log_mdd=fix_point_strategy(lmdd0);
    
else
    
    log_mdd=iterative_strategy(lmdd0);
    
end

extras=struct('bridge_lmdd',lmdd,'bridge_conv',conv);

    function [log_mdd]=iterative_strategy(lmdd0)
        % 3- for a number of iterations, update an initial guess for the
        % log MDD
        %----------------------------------------------------------------
        
        lmdd(1)=lmdd0;
        
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
            
            if opts.debug
                
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

    function log_mdd=fix_point_strategy(lmdd0)
        
        fzero_options=optimset('Display','none');
        
        if opts.debug
            
            fprintf('BRIDGE: Now bisecting...');
            
            tic
            
            fzero_options=optimset('Display','iter');
            
        end
        
        next=0;
        
        log_mdd=fzero(@one_right_hand_side,lmdd0,fzero_options);
        
        if opts.debug
            
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
            exp(log(opts.L)+f1)+...
            exp(log(obj.M)+f2)...
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

    function set_inputs()
        
        if n<3
            
            opts=[];
            
            if n<2
                
                fix_point=[];
                
            end
            
        end
        
        if isempty(fix_point),fix_point=true; end
        
        assert(islogical(fix_point) && isscalar(fix_point),...
            'fix_point must be true or alse')
        
        opts=mdd.global_options(opts);
        
    end

end
