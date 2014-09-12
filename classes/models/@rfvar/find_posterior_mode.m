function [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb,ub)

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
    
    [vdata,bigx,bigy]=restricted_var_data(obj);
    a_prior=[obj.estimation.priors.prior_mean];
    a_prior=a_prior(vdata.estim_locs);
    a_prior=a_prior(:);
    va_prior=[obj.estimation.priors.prior_stdev];
    va_prior=diag(va_prior(vdata.estim_locs).^2);
    
    a2tilde_prior=struct('a',obj.linear_restrictions_data.a2tilde_func(a_prior),...
        'V',obj.linear_restrictions_data.a2tilde_func(va_prior,true));
    
    [a2tilde_post,a2tilde_ols]=posterior_mode_engine(vdata.Xtilde,...
        vdata.ytilde,a2tilde_prior);
    
    resid_post=reshape(a2tilde_post.resid,obj.endogenous.number(end),[]);
    [nvars,nobs]=size(resid_post);
    K=(npar-sum(1:nvars))/nvars;
    vcov=(resid_post*resid_post.')/(nobs-K);
    
    resid_ols=reshape(a2tilde_ols.resid,obj.endogenous.number(end),[]);
    a2tilde_ols.SSE=(resid_ols*resid_ols.');
    a2tilde_ols.SIGMA = a2tilde_ols.SSE./(nobs-K);
    
    a_post=obj.linear_restrictions_data.a_func(a2tilde_post.a);
    
    x1=vartools.build_parameter_vector(vdata,a_post,vcov);
    
    a2tilde.prior=a2tilde_prior;
    a2tilde.ols=a2tilde_ols;
    a2tilde.post=a2tilde_post;
    
    f1=uminus(log_posterior_kernel(obj,x1));
    f0=uminus(log_posterior_kernel(obj,x0));
    warning([mfilename,':: this wrong variance for one-regime VARs has to change'])
    H=eye(npar);
    viol=[];
    funevals=1;
    issue={};
    
    % prepare for subsequent simulation
    %----------------------------------
    
    a2Aprime=@(x)transpose(reshape(obj.linear_restrictions_data.a_func(x),nvars,K));
    
    if any(strcmp(obj.options.vp_prior_type,{'normal_wishart','indep_normal_wishart'}))
        % Hyperparameters on inv(SIGMA) ~ W(prior.dof_SIGMA,inv(prior.scale_SIGMA))
        a2tilde.prior.dof_SIGMA = endo_nbr+1;         %<---- prior Degrees of Freedom (DoF) of SIGMA
        a2tilde.prior.scale_SIGMA = eye(endo_nbr);    %<---- prior scale of SIGMA
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
        a2tilde.post.dof_SIGMA = nobs + a2tilde.prior.dof_SIGMA;
        if strcmp(prior_type,'normal_wishart')
            % we invert and then apply the function:probably the simplest thing
            % to do
            %------------------------------------------------------------------
            iVa=a_func(inv(a2tilde.prior.V),true);
            XpX=bigx*bigx';
            A_OLS=a2Aprime(a2tilde.ols.a);
            a2tilde.post.scale_SIGMA = a2tilde.ols.SSE + a2tilde.prior.scale_SIGMA + ...
                A_OLS'*XpX*A_OLS + ...
                A_prior'*iVa*A_prior - ...
                A_post'*(iVa + XpX)*A_post;
        end
    end
    
    obj.constant_var_data=struct('a2tilde',a2tilde,...
        'a_func',obj.linear_restrictions_data.a_func,...
        'K',K,...
        'endo_nbr',nvars,...
        'nobs',nobs,...
        'X',transpose(bigx),...
        'Y',transpose(bigy),...
        'prior_type',obj.options.vp_prior_type,...
        'SIGMA',a2tilde_ols.SIGMA,...
        'x0',x1,...
        'f0',f1,...
        'funevals',0,...
        'a2Aprime',a2Aprime,...
        'na2',numel(a2tilde.prior.a),...
        'vdata',vdata);
else
    [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode@rise_generic(obj,x0,lb,ub);
end

    function [post,ols]=posterior_mode_engine(X,y,prior)
        
        npar_short=obj.linear_restrictions_data.npar_short;
        iVa_prior=prior.V\eye(npar_short);
        
        % ols
        %----
        ols.a=X\y;
        ols.resid=y-X*ols.a;
        iV_ols=(X'*X);
        
        % posterior
        %----------
        post.V=(iV_ols+iVa_prior)\eye(npar_short);
        post.a=post.V*(iV_ols*ols.a+iVa_prior*prior.a);
        post.resid=y-X*post.a;
    end
end

function [vd,bigx,bigy]=restricted_var_data(obj)
[bigy,bigx,nv]=vartools.set_y_and_x(obj.data.y,obj.data.x,...
    obj.nlags,obj.constant);
vd=struct();
xi=kron(bigx',eye(nv));
f=obj.linear_restrictions_data.R1i_r_0;
G=obj.linear_restrictions_data.R1i_R2_I2;
vd.ytilde=bigy(:);
if any(f)
    vd.ytilde=vd.ytilde-xi*f(obj.linear_restrictions_data.ievec);
end
vd.Xtilde=xi*G(obj.linear_restrictions_data.ievec,:);

vd.orig_estim_names={obj.estimation.priors.name};
% vd.estim_names=estim_names;
vd.estim_locs=locate_variables(obj.linear_restrictions_data.estim_names,vd.orig_estim_names);
end
%}
