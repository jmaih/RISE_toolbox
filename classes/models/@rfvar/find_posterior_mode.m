function [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb,ub)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


nobj=numel(obj);
if nobj==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    x1=struct();
    return
end

if obj.markov_chains.regimes_number==1 && obj.options.vp_analytical_post_mode
    % the following is hard-coded
    %----------------------------
    compute_hessian=false;
    
    npar=size(x0,1);
    a2tilde=struct();
    
    [vdata,bigx,bigy,orig_order]=restricted_var_data(obj);
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
    
    fh=@(x)uminus(log_posterior_kernel(obj,x));
    f1=fh(x1);
    f0=fh(x0);
    % compute Hessian
    %----------------
    if compute_hessian
        switch lower(obj.options.hessian_type)
            case 'fd'
                H = utils.hessian.finite_differences(fh,x1);
            case 'opg'
                H = utils.hessian.outer_product(fh,x1);
                if any(any(isnan(H)))||any(any(isinf(H)))
                    issue='OPG unstable and inaccurate for calculation of Hessian, switched to finite differences';
                    warning([mfilename,':: ',issue]) %#ok<WNTAG>
                    warning([mfilename,':: OPG unstable for calculation of Hessian, switching to finite differences']) %#ok<WNTAG>
                    H = finite_difference_hessian(fh,x1);
                end
            otherwise
                issue=['unknow hessian option ',hessian_type,' using finite differences'];
                warning([mfilename,':: ',issue]) %#ok<WNTAG>
                H = finite_difference_hessian(fh,x1);
        end
    else
        warning([mfilename,':: this variance for one-regime VARs is wrong'])
        H=eye(npar);
    end
    % add remaining items
    %--------------------
    viol=[];
    funevals=2;
    issue='';%issue={};
    
    % prepare for subsequent simulation
    %----------------------------------
    
    a2Aprime=@(x)transform_to_matrix_form(obj.linear_restrictions_data.a_func(x));
    
    if any(strcmp(obj.options.vp_prior_type,{'normal_wishart','indep_normal_wishart'}))
        % Hyperparameters on inv(SIGMA) ~ W(prior.dof_SIGMA,inv(prior.scale_SIGMA))
        a2tilde.prior.dof_SIGMA = nvars+1;         %<---- prior Degrees of Freedom (DoF) of SIGMA
        a2tilde.prior.scale_SIGMA = eye(nvars);    %<---- prior scale of SIGMA
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
        a2tilde.post.dof_SIGMA = nobs + a2tilde.prior.dof_SIGMA;
        if strcmp(obj.options.vp_prior_type,'normal_wishart')
            % we invert and then apply the function:probably the simplest thing
            % to do
            %------------------------------------------------------------------
            iVa=obj.linear_restrictions_data.a_func(inv(a2tilde.prior.V),true);
            XpX=bigx*bigx';
            A_OLS=a2Aprime(a2tilde.ols.a);
            A_prior=a2Aprime(a2tilde.prior.a);
            A_post=a2Aprime(a2tilde.post.a);
            % %             iVa=reshape(diag(iVa),size(A_post));
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

    function A=transform_to_matrix_form(x)
        x=x(orig_order);
        A=reshape(x(:),nvars,K);
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

function [vd,bigx,bigy,orig_order]=restricted_var_data(obj)
[bigy,bigx,nv]=vartools.set_y_and_x(obj.data.y,obj.data.x,...
    obj.nlags,obj.constant);
vd=struct();
xi=kron(bigx',eye(nv));
% re-order the columns of xi to conform with the order of the names of the
% estimated parameters
inv_order=locate_variables(obj.all_param_names_vec,obj.linear_restrictions_data.estim_names,true);
good=~isnan(inv_order);
% tmp=obj.all_param_names_vec(good);
inv_order=inv_order(good);
if numel(inv_order)~=numel(obj.linear_restrictions_data.estim_names)
    error('Please contact junior.maih@gmail.com with this')
end
% this is expected to be perfectly symmetric such that orig_order =
% inv_order
orig_order(inv_order)=1:numel(inv_order);
xi=xi(:,inv_order);

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
vd.not_estimated=setdiff(1:numel(vd.orig_estim_names),vd.estim_locs);

pp=regexp(vd.orig_estim_names(vd.not_estimated),...
    '\w+(?<first>\d+)_(?<second>\d+)','names');
pp=[pp{:}];
vd.p1=cell2mat(cellfun(@(x)str2double(x),{pp.first},'uniformOutput',false));
vd.p2=cell2mat(cellfun(@(x)str2double(x),{pp.second},'uniformOutput',false));
vd.same=vd.p1==vd.p2;
end
%}
