function [LogLik,Incr,retcode,obj]=likelihood_smooth_transition_var(params,obj)

% parameters are as follows
% p*(nlags*p+1+nx)*nregs unrestricted
% p standard deviations >0
% p*(p-1)/2 correlations [-1,1]
% p*(nregs-1) g parameters >= 0
% threshold parameters that must lie in the boundaries of the threshold
% variable and sorted...

LogLik=-1e+8;

Incr=[];

obj=assign_estimates(obj,params); 

[u,~,retcode,obj]=low_level_residuals(obj);

if retcode
    
    return
    
end

s=obj.solution;

SIG=s.C*s.C.';

det_SIG=det(SIG);

is_failed=det_SIG<=1e-10;% ~all(isfinite(iSIG(:)));

if is_failed
    
    retcode=30002;
    
    return
    
end

p=obj.endogenous.number;

ldet_pl2pi=log(det_SIG)+p*log(2*pi);

iSIG=SIG\eye(p);

Incr=-.5*(diag((u.'*iSIG*u))+ldet_pl2pi);

LogLik=sum(Incr);

if ~utils.error.valid(LogLik)
    
    % nans or inf in likelihood
    retcode=300002;
    
end

end
