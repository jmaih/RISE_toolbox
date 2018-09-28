function [log_mdd,extras]=mueller(obj,~,opts)

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

% Logq_L: weighting draws evaluated through weighting function
%--------------------------------------------------------------

% LogPost_L: weighting draws evaluated through the log-posterior kernel
%----------------------------------------------------------------
[LogPost_L,Logq_L]=iid_draws(obj,iid_draws_,opts,'MUELLER');

% old draws evaluated through the weighting function
%----------------------------------------------------
Logq_M=old_draws_in_weighting_function(obj,opts,'MUELLER');

Logq_M=Logq_M(1:nold);

LogPost_M=obj.LogPost_M(1:nold);

% get the ratios
%----------------
log_rh=LogPost_L-Logq_L;

log_rp=Logq_M-LogPost_M;

fzero_options=optimset('Display','none');

if opts.debug
    
    fprintf('MUELLER: Now bisecting...');
    
    tic
    
    fzero_options=optimset('Display','iter');
    
end

next=0;

max_iter=1000;

record=zeros(max_iter,1);

log_mdd=fzero(@hfunc,0,fzero_options);

if opts.debug
    
    fprintf('Done in %0.4d seconds\n',toc);
    
end

extras=struct('mueller_convergence',record(1:next));

log_mdd=obj.thecoef*log_mdd;

    function h=hfunc(log_c)
        % left side
        %------------
        left=log_c+log_rh;
        
        good=left<0;
        
        left=sum(1-exp(left(good)))/opts.L;
        
        % right side
        %------------
        right=log_rp-log_c;
        
        good=right<0;
        
        right=sum(1-exp(right(good)))/nold;
        
        h=left-right;
        
        next=next+1;
        
        if next==max_iter
            
            record(next+1000)=0;
            
            max_iter=max_iter+1000;
            
        end
        
        record(next)=h;
        
    end

    function set_inputs()
        
        if n<3
            
            opts=[];
            
        end
        
        opts=mdd.global_options(opts);
        
    end

end
