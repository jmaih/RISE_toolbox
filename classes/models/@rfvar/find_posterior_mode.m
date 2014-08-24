function [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb,ub,...
    basics,general_restrictions,gen_restr_args)

nobj=numel(obj);
if nobj==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    x1=struct();
    return
end

if obj.markov_chains.regimes_number==1 && obj.options.vp_analytical_post_mode
    npar=size(x0,1);
    a2tilde=struct();
    
    vdata=restricted_var_data(obj,basics);
    a_prior=[obj.estimation.priors.prior_mean];
    a_prior=a_prior(vdata.estim_locs);
    a_prior=a_prior(:);
    va_prior=[obj.estimation.priors.prior_stdev];
    va_prior=diag(va_prior(vdata.estim_locs).^2);
    
    a2tilde_prior=struct('a',basics.a2tilde_func(a_prior),...
        'V',basics.a2tilde_func(va_prior,true));
    
    [a2tilde_post,a2tilde_ols]=posterior_mode_engine(...
        vdata.Xtilde,vdata.ytilde,a2tilde_prior);
    
    resid_post=reshape(a2tilde_post.resid,obj.endogenous.number(end),[]);
    [nvars,nobs]=size(resid_post);
    K=(npar-sum(1:nvars))/nvars;
    vcov=(resid_post*resid_post.')/(nobs-K);
    
    resid_ols=reshape(a2tilde_ols.resid,obj.endogenous.number(end),[]);
    a2tilde_ols.SSE=(resid_ols*resid_ols.');
    a2tilde_ols.SIGMA = a2tilde_ols.SSE./(nobs-K);
    
    a_post=basics.a_func(a2tilde_post.a);
    
    x1=vartools.build_parameter_vector(vdata,a_post,vcov);
    
    warning([mfilename,':: this wrong variance for one-regime VARs has to change'])
    a2tilde.prior=a2tilde_prior;
    a2tilde.ols=a2tilde_ols;
    a2tilde.post=a2tilde_post;
    
    f1=uminus(log_posterior_kernel(obj,x1));
    f0=uminus(log_posterior_kernel(obj,x0));
    H=eye(npar);
    viol=[];
    funevals=1;
    issue={};
    keyboard
else
    [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode@rise_generic(obj,x0,lb,ub,...
    basics,general_restrictions,gen_restr_args);
end

    function [post,ols]=posterior_mode_engine(X,y,prior)
        
        iVa_prior=prior.V\eye(basics.npar_short);
        
        % ols
        %----
        ols.a=X\y;
        ols.resid=y-X*ols.a;
        iV_ols=(X'*X);
        
        % posterior
        %----------
        iV_post=(iV_ols+iVa_prior);
        post.a=iV_post\(iV_ols*ols.a+iVa_prior*prior.a);
        post.resid=y-X*post.a;
    end

end

function vd=restricted_var_data(obj,basics)
[bigy,bigx,nv]=vartools.set_y_and_x(obj.data.y,obj.data.x,...
    obj.nlags,obj.constant);
vd=struct();
xi=kron(bigx',eye(nv));
f=basics.R1i_r_0;
G=basics.R1i_R2_I2;
vd.ytilde=bigy(:);
if any(f)
    vd.ytilde=vd.ytilde-xi*f(basics.ievec);
end
vd.Xtilde=xi*G(basics.ievec,:);

vd.orig_estim_names={obj.estimation.priors.name};
% vd.estim_names=estim_names;
vd.estim_locs=locate_variables(basics.estim_names,vd.orig_estim_names);
end
%}
