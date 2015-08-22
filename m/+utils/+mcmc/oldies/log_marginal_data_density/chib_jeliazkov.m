classdef chib_jeliazkov
    % chib_jeliazkov computation of the MDD using chib and jeliazkov
    % - [chib_jeliazkov](chib_jeliazkov/chib_jeliazkov)
    % - [conclude](chib_jeliazkov/conclude)
    % - [update](chib_jeliazkov/update)
    properties
        xmode
        vcov_mode
        log_post_mode
        lower_bound
        upper_bound
        log_post_func
    end
    properties(Access=protected)
        det_vcov
        vcov_inv
        npar
        CJ=0;
        sample_size=0;
    end
    methods
        function obj=chib_jeliazkov(xmode,vcov_mode,log_post_mode,...
                lower_bound,upper_bound,log_post_func)
            obj.xmode=xmode;
            obj.npar=numel(xmode);
            obj.det_vcov=det(vcov_mode);
            obj.vcov_mode=vcov_mode;
            obj.vcov_inv=vcov_mode\eye(obj.npar);
            obj.lower_bound=lower_bound;
            obj.upper_bound=upper_bound;
            obj.log_post_mode=log_post_mode;
            obj.log_post_func=log_post_func;
        end
        function obj=update(obj,param,~)
            lik=conditional_likelihood(obj.xmode-param,obj.vcov_inv,obj.det_vcov,obj.npar);
            obj.CJ=obj.CJ+lik;
            obj.sample_size=obj.sample_size+1;
        end
        function log_mdd=conclude(obj,varargin)
            obj.CJ=obj.CJ/obj.sample_size;
            J=obj.sample_size;
            CS=transpose(chol(obj.vcov_mode));
            f0=obj.log_post_mode;
            alpha0=0;
            x_mode=obj.xmode;
            for ii=1:J
                valid=false;
                while ~valid
                    xi=x_mode+CS*randn(obj.npar,1);
                    if all(xi>=obj.lower_bound) && all(xi<=obj.upper_bound)
                        valid=true;
                    end
                end
                f_i=obj.log_post_func(xi);
                a=utils.mcmc.alpha_probability(f_i,f0);
                alpha0=alpha0+a;
            end
            alpha0=alpha0/J;
            log_mdd=log(obj.CJ/alpha0);
        end
    end
end
