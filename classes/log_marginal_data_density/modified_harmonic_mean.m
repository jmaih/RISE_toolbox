classdef modified_harmonic_mean
    % modified_harmonic_mean computation of the MDD using modified harmonic mean
    % - [conclude](modified_harmonic_mean/conclude)
    % - [modified_harmonic_mean](modified_harmonic_mean/modified_harmonic_mean)
    % - [update](modified_harmonic_mean/update)
    properties
        xmean
        vcov_mean
    end
    properties(Access=protected)
        det_vcov
        vcov_inv
        npar
        MHM=0;
        sample_size=0;
    end
    methods
        function obj=modified_harmonic_mean(xmean,vcov_mean,~,~,~,~)
            obj.xmean=xmean;
            obj.npar=numel(xmean);
            obj.det_vcov=det(vcov_mean);
            obj.vcov_mean=vcov_mean;
            obj.vcov_inv=vcov_mean\eye(obj.npar);
        end
        function obj=update(obj,param,log_post_param)
            lnf=geweke_truncated_multivariate_normal(param,obj.xmean,obj.det_vcov,obj.vcov_inv);
            obj.MHM=obj.MHM+exp(lnf-log_post_param);
            obj.sample_size=obj.sample_size+1;
        end
        function log_mdd=conclude(obj,varargin)
            obj.MHM=1./(obj.MHM/obj.sample_size);
            log_mdd=max(log(obj.MHM));
        end
    end
end

function lnf=geweke_truncated_multivariate_normal(theta,theta_mean,detV,Vi,tau)
if nargin<5
    tau=[];
end
if isempty(tau)
    tau=transpose(0.1:.1:.9);
end
d=numel(theta);
[lik,v_iF_v]=conditional_likelihood(theta-theta_mean,Vi,detV,d);
f=lik./tau.*(v_iF_v<=chi2inv(tau,d));
lnf=log(f);
end