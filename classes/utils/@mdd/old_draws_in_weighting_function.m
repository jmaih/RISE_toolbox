function Logq_M=old_draws_in_weighting_function(obj,opts,caller)

if opts.debug
    
    fprintf('%s: Now evaluating %0.0d old draws in weighting function...',caller,obj.M);
    
    tic
    
end

moments(obj,opts)

m=obj.moms.m;

nw=@(v)mdd.normal_weighting(v,m.det_Sigma,m.Sigma_i);

Logq_M=zeros(1,obj.M);

for idraw=1:obj.M
    
    draw=obj.theta_draws(:,idraw);
    
    v=draw-m.x_bar;
    
    Logq_M(idraw)=nw(v);
    
end

if opts.debug
    
    fprintf('Done in %0.4d seconds\n\n',toc);
    
end

end