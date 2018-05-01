function [LogLik,Incr,retcode,obj]=likelihood_dsge_var(params,obj)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% N.B: This function assumes the likelihood is to be maximized

if nargin~=2
    
    error([mfilename,':: Number of arguments must be 0 or 2'])
    
end

% this important output is never created
Incr=[];

obj=assign_estimates(obj,params);

[obj,retcode]=solve(obj);

LogLik=-obj.options.estim_penalty;

if ~retcode
    
    p=obj.dsge_var.p;% the var order
    
    T=obj.dsge_var.T;% the sample size
    
    const=obj.dsge_var.constant; % flag for constant
    
    n=obj.dsge_var.n; % number of variables
    
    k=const+n*p;
    
    SIG_theta=obj.dsge_var.var_approx.SIG;
    
    GXX=obj.dsge_var.var_approx.GXX;
    
    YY=obj.dsge_var.YY;
    
    YX=obj.dsge_var.YX;
    
    XX=obj.dsge_var.XX;
    
    XY=obj.dsge_var.XY;
    
    PHIb=obj.dsge_var.posterior.PHI;
    
    SIGb=obj.dsge_var.posterior.SIG;
    
    ltgxx=obj.dsge_var.posterior.ZZ;
    
    lambda=obj.dsge_var.lambda;
    
    if isinf(lambda)
        
        LogLik = -.5*T*log(det(SIGb))...
            -.5*n*T*log(2*pi)...
            -.5*trace(SIGb\(YY-YX*PHIb-PHIb'*XY+PHIb'*XX*PHIb));
        
    else
        
        LogLik=-.5*n*log(det(ltgxx))...
            -.5*((1+lambda)*T-k)*log(det((1+lambda)*T*SIGb))...
            +.5*n*log(det(lambda*T*GXX))...
            +.5*(lambda*T-k)*log(det(lambda*T*SIG_theta))...
            -.5*n*T*log(2*pi)...
            +.5*n*T*log(2)...
            +sum(gammaln(.5*((1+lambda)*T-k+1-(1:n))))...
            -sum(gammaln(.5*(lambda*T-k+1-(1:n))));
        
    end
    
    if isnan(LogLik)||isinf(LogLik)||~isreal(LogLik)
        
        retcode=300002;
        
    end
    
    if obj.options.kf_filtering_level
        % now we filter the data, provided, the parameters
        % estimated using the dsge-var do not have a low density as
        % from the point of view of the pure dsge.
        [obj,dsge_log_lik,~,dsge_retcode]=filter(obj);
        
        LogLik=[LogLik,dsge_log_lik];
        
        retcode=[retcode,dsge_retcode];
        
    end
    
end

if obj.options.debug
    
    disp(LogLik)
    
end

end

