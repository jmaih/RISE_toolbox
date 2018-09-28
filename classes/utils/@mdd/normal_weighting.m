function [loglik,v_iF_v,lik]=normal_weighting(v,det_Sigma,Sigma_i)

v_iF_v=v'*Sigma_i*v;

d=size(v,1);

loglik=-.5*(d*log(2*pi)+log(det_Sigma)+v_iF_v);

if nargout>2
    
    lik=exp(loglik);
    
end

end
