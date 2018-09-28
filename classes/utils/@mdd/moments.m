function moments(obj,opts)

if isempty(opts.center_at_mean)
    
    opts.center_at_mean=false;
    
end

if opts.debug
    
    fprintf('MDD: Now computing moments from %0.0d draws...\n',obj.M);
    
    tic
    
end

if ~isempty(obj.moms) && obj.moms.center_at_mean~=opts.center_at_mean
    
    return
    
end

m=struct();

if opts.center_at_mean
    
    m.x_bar = (1/obj.M)*sum(obj.theta_draws,2);
    
else
    
    m.x_bar = obj.theta_mode.x;
    
end

S=obj.theta_draws-m.x_bar(:,ones(1,obj.M));

m.Sigma_bar=S*S.'/obj.M; %<---m.Sigma_bar =cov(obj.theta_draws.');

m.Sigma_i=m.Sigma_bar\eye(obj.d);

m.det_Sigma=det(m.Sigma_bar);

m.Shat = chol(m.Sigma_bar,'lower');

if opts.debug
    
    fprintf('Done in %0.4d seconds\n\n',toc);
    
end

obj.moms.m=m;

obj.moms.center_at_mean=opts.center_at_mean;

end
