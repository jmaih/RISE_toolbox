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
% More About
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

if obj.markov_chains.regimes_number==1 && obj.options.estim_analytical_post_mode
    % the following is hard-coded
    %----------------------------
    compute_hessian=false;
    
    npar=size(x0,1);
    a2tilde=struct();
    
    [vdata,bigx,bigy,orig_order,inv_order,smpl]=restricted_var_data(obj);
    [nvars,nobs]=size(bigy);
    a_prior=[obj.estimation.priors.prior_mean];
    a_prior=a_prior(vdata.estim_locs);
    a_prior=a_prior(:);
    va_prior=[obj.estimation.priors.prior_stdev];
    va_prior=diag(va_prior(inv_order).^2);% <---va_prior=diag(va_prior(vdata.estim_locs).^2);
    
    a2tilde_prior=struct('a',obj.linear_restrictions_data.a2tilde_func(a_prior),...
        'V',obj.linear_restrictions_data.a2tilde_func(va_prior,true),...
        'dof_SIGMA',nvars+1,...%<---- prior Degrees of Freedom (DoF) of SIGMA
        'scale_SIGMA',eye(nvars)...%<---- prior scale of SIGMA
        );
        % Hyperparameters on inv(SIGMA) ~ W(prior.dof_SIGMA,inv(prior.scale_SIGMA))
    
    [a2tilde_post,a2tilde_ols,estimafy]=posterior_mode_engine(vdata.Xtilde,...
        vdata.ytilde,a2tilde_prior);
    
    resid_post=reshape(a2tilde_post.resid,obj.endogenous.number,[]);
    K=(npar-sum(1:nvars))/nvars;
    vcov=(resid_post*resid_post.')/(nobs-K);
    
    resid_ols=reshape(a2tilde_ols.resid,obj.endogenous.number,[]);
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
    % compute Hessian: we always compute two hessians, one returned by the optimizer and the other computed numerically
    %----------------
    H=repmat(eye(npar,npar),[1,1,2]);
    if compute_hessian
		% compute the numerical hessian
		%-------------------------------
		[H(:,:,2),issue]=utils.hessian.numerical(fh,x1,lower(obj.options.hessian_type));
    else
        warning([mfilename,':: this variance for one-regime VARs is wrong'])
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
        if strcmp(obj.options.vp_prior_type,'normal_wishart')
            % we invert and then apply the function:probably the simplest thing
            % to do
            %------------------------------------------------------------------
            % unrestrict and then take the first K guys
            iVa=diag(1./diag(va_prior(1:K,1:K)));
            XpX=bigx*bigx';
            A_OLS=a2Aprime(a2tilde.ols.a);
            A_prior=a2Aprime(a2tilde.prior.a);
            A_post=a2Aprime(a2tilde.post.a);
            a2tilde.post.scale_SIGMA = a2tilde.ols.SSE + a2tilde.prior.scale_SIGMA + ...
                A_OLS*XpX*A_OLS' + ...
                A_prior*iVa*A_prior' - ...
                A_post*(iVa + XpX)*A_post';
        end
    end
    
    obj.constant_var_data=struct('a2tilde',a2tilde,...
        'a_func',obj.linear_restrictions_data.a_func,...
        'a2tilde_func',obj.linear_restrictions_data.a2tilde_func,...
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
        'vdata',vdata,...
        'estimafy',estimafy,...
        'inv_order',inv_order);
else
    [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode@rise_generic(obj,x0,lb,ub);
end

    function A=transform_to_matrix_form(x)
        x=x(orig_order);
        A=reshape(x(:),nvars,K);
    end

    function [post,ols,estimafy]=posterior_mode_engine(X,y,prior)
        
        npar_short=obj.linear_restrictions_data.npar_short;
        iVa_prior=prior.V\eye(npar_short);
        
        % gls
        %----
        iSIGU=1;
        switch obj.options.vp_prior_type
            case {'minnesota','indep_normal_wishart'}
                if obj.options.vp_gls_ar1_processes
                    iSIGU=diag(1./diag(obj.miscellaneous.constant_var.sigma_));
                else
                    results=vartools.ols(bigy,bigx,0,false);
                    iSIGU=results.SIGols\eye(obj.endogenous.number);%
                end
                iSIGU=kron(eye(smpl),iSIGU);%<---iSIGU=kron(iSIGU,eye(smpl));
                % N.B: for the indep_normal_wishart, This is not the exact
                % formula but we still need to start somewhere for the
                % initialization of the Gibbs sampler.
            case 'none'
                iVa_prior(:)=0;
            case 'normal_wishart'
            case {'jeffrey','diffuse'}
            otherwise
        end
        
        % posterior and ols
        %-------------------
        [post,ols]=compute_posterior_and_ols(iSIGU);
        
        estimafy=@compute_posterior_and_ols;
        
        function [post,ols]=compute_posterior_and_ols(iSIGU,flag)
            if nargin<2
                flag=false;
            end
            if flag
                iSIGU=kron(eye(smpl),iSIGU);
            end
            iV_ols=(X'*iSIGU*X);
            ols.a=(X'*iSIGU*X)\(X'*iSIGU*y);
            ols.resid=y-X*ols.a;
            post.V=(iV_ols+iVa_prior)\eye(npar_short);
            post.a=post.V*(iV_ols*ols.a+iVa_prior*prior.a);
            post.resid=y-X*post.a;
            % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
            post.dof_SIGMA = nobs + prior.dof_SIGMA;
% %             % fishy business : we use different variances for the computation
% %             % of the mean and for posterior simulation
% %             %------------------------------------------------------------------
% %             if strcmp(obj.options.vp_prior_type,'jeffrey')
% %                 results=vartools.ols(bigy,bigx,0,false);
% %                 iSIGU_fishy=results.SIGols\eye(obj.endogenous.number);%
% %                 post.V=inv(X'*kron(eye(smpl),iSIGU_fishy)*X);
% %             end
        end
    end

end

function [vd,bigx,bigy,orig_order,inv_order,smpl]=restricted_var_data(obj)
[bigy,bigx,nv,smpl]=vartools.set_y_and_x(obj.data.y,obj.data.x,...
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
