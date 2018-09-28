function [log_mdd,extras]=cj(obj,~,opts)

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

% Posterior kernel at the mode
%--------------------------------------------------
numerator=max(obj.LogPost_M);

% weighting-function density for the M original draws
%-----------------------------------------------------
Logq_M=old_draws_in_weighting_function(obj,opts,'CJ');

Logq_M=Logq_M(1:nold);

% log posterior kernel of the IID draws
%---------------------------------------
[LogPost_L]=iid_draws(obj,iid_draws_,opts,'CJ');

% acceptance probability for the IID draws relative to the
% posterior mode
%-----------------------------------------------------------------
alpha=zeros(1,opts.L);

for idraw=1:opts.L
    
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

    function set_inputs()
        
        if n<3
            
            opts=[];
            
        end
        
        opts=mdd.global_options(opts);
        
    end

end
