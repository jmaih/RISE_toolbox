function log_mhm=mhm(obj,mhm_tau,opts)

n=nargin;

set_options()

critical_value=chi2inv(mhm_tau,obj.d);

log_tau=log(mhm_tau);

post_tmn=zeros(1,numel(log_tau));

if opts.debug
    
    fprintf('MHM: Just doing it...');
    
    tic
    
end

moments(obj,opts)

m=obj.moms.m;

nw=@(v)mdd.normal_weighting(v,m.det_Sigma,m.Sigma_i);

% we serialize the following in order to save on memory and because
% calling conditional_likelihood should be relatively cheap
for idraw=1:obj.M % parfor (idraw=1:obj.M,NumWorkers)
    % evaluate proposal
    %------------------
    [loglik,v_iF_v]=nw(obj.theta_draws(:,idraw)-m.x_bar);
    
    % geweke's truncated multivariate normal
    %---------------------------------------
    good=v_iF_v<=critical_value;
    
    log_f=loglik-log_tau(good);
    
    post_tmn(good)=post_tmn(good)+exp(log_f-obj.LogPost_M(idraw)+obj.facilitator);
    
end

if opts.debug
    
    fprintf('Done in %0.4d seconds\n',toc);
    
end

post_tmn=post_tmn/obj.M;

log_mhm=max(obj.facilitator-log(post_tmn));

    function set_options()
        
        if n<3
            
            opts=[];
            
            if nargin<2
                
                mhm_tau=[];
                
            end
            
        end
        
        if isempty(mhm_tau)
            
            mhm_tau=(.1:.1:.9);
            
        end
        
        assert(all(mhm_tau>0) && all(mhm_tau<1),...
            'mhm_tau (truncation probabilities of MHM) must all be in (0,1)')
        
        opts=mdd.global_options(opts);
        
    end

end
