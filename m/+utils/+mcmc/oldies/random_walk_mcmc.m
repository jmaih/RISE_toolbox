function [x1,f1,accepted,funevals,alpha_prob]=random_walk_mcmc(...
minus_log_post_func,x0,f0,drawfun,cCS,mcmc_delay_rejection,funevals)
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

[theta_s,minusLogPost_s]=new_proposal();
alpha_prob=utils.mcmc.alpha_probability(-minusLogPost_s,-f0);
accepted=alpha_prob>rand;
if accepted
    x1=theta_s;
    f1=minusLogPost_s;
else
    x1=x0;
    f1=f0;
    if mcmc_delay_rejection
        [theta_s2,minusLogPost_s2]=new_proposal();
        alpha_prob2=utils.mcmc.alpha_probability(-minusLogPost_s2,-minusLogPost_s);
        alpha13 = exp(-minusLogPost_s2-(-f0))*...
            qdens(theta_s2,theta_s)/qdens(x0,theta_s)*...
            (1-alpha_prob2)/(1-alpha_prob);
        accepted=alpha13>rand;
        if accepted
            x1=theta_s2;
            f1=minusLogPost_s2;
        end
    end
end
    function qab=qdens(a,b)
        ab=a-b;
        ab=ab(:);
        qab=exp(-0.5*(ab'*ab));
    end
    function [d,minusLogPost]=new_proposal()
        d=drawfun(x0,cCS);% d=x0+cCS*randn(npar,1);
        minusLogPost=minus_log_post_func(d);
        funevals=funevals+1;
    end
end