function [LogPost_L,Logq_L]=iid_draws(obj,iid_draws_,opts,caller)

if opts.debug
    
    fprintf('%s: Now evaluating %0.0d IID draws...',caller,opts.L);
    
    tic
    
end

eval_weighting=nargout>1;

if eval_weighting
    
    Logq_L=zeros(1,opts.L);
    
end

is_drawn=~isempty(iid_draws_);

moments(obj,opts)

m=obj.moms.m;

nw=@(v)mdd.normal_weighting(v,m.det_Sigma,m.Sigma_i);

x_bar=m.x_bar;

Shat=m.Shat;

d=numel(x_bar);

for idraw=1:opts.L
    
    if ~is_drawn
        
        if idraw==1
            
            iid_draws_=struct();
            
        end
        
        v=Shat*randn(d,1);
        
        draw=obj.recenter(x_bar+v);
        
        iid_draws_(idraw).x=draw;
        
        iid_draws_(idraw).f=obj.log_post_kern(draw);
        
    end
    
    if eval_weighting
        
        v=iid_draws_(idraw).x-x_bar;
        
        Logq_L(idraw)=nw(v);
        
    end
    
end

LogPost_L=obj.thecoef*[iid_draws_.f];

if opts.debug
    
    fprintf('Done in %0.4d seconds\n\n',toc);
    
end

end
