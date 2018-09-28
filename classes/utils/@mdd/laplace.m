function log_mdd=laplace(obj,~,~)

[lpbest,loc]=max(obj.LogPost_M);

theta_mode=obj.theta_draws(:,loc);

H=utils.hessian.finite_differences(obj.log_post_kern,theta_mode);

if obj.maximization
    
    H=-H;
    
end

Hinv=H\eye(size(H,1));

log_mdd=utils.marginal_data_density.laplace_mdd(lpbest,Hinv);

end
