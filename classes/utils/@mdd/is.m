function log_mdd=is(obj,~,opts)

% Sylvia Frühwirth-Schnatter (2004): "Estimating marginal
% likelihoods for mixture and Markov switching models using bridge
% sampling techniques". Econometrics Journal 7(1) pp 143--167

n=nargin;

set_inputs()

iid_draws_=[];

c=obj.thecoef;

if opts.draws_are_iid
    
    iid_draws_=struct();
    
    for idraw=1:opts.L
        
        iid_draws_(idraw).x=obj.theta_draws(:,end-idraw+1);
        
        % If maximization, c=1. If minimization, c=-1 and the current
        % operation will be reversed in iid_draws;
        iid_draws_(idraw).f=c*obj.LogPost_M(end-idraw+1);
        
    end
    
end

[LogPost_L,Logq_L]=iid_draws(obj,iid_draws_,opts,'IS');

ratio = LogPost_L-Logq_L;

ratiomax=max(ratio);

log_mdd = ratiomax+log(mean(exp(ratio-ratiomax)));

    function set_inputs()
        
        if n<3
            
            opts=[];
            
        end
        
        opts=mdd.global_options(opts);
        
    end

end
