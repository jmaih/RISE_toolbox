function [smpl,fsmpl,accept_rate,start] = constant_bvar_sampler(start,nsamples,varargin)
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


pnames_defaults = {
    'burnin',0,@(x)isa(x,'double') && isscalar(x) && x==round(x) && x>=0
    'thin',1,@(x)isa(x,'double') && isscalar(x) && x==round(x) && x>=1
    'waitbar_update',[],@(x)isa(x,'function_handle')
    };
if nargin==0
    smpl=pnames_defaults(:,1:3);
    return
end

%---------------------------------------------
prior_type=start.prior_type;
Y=start.Y ;
X=start.X ;
nobs=start.nobs ;
K=start.K ;
a_func=start.a_func ;
a2tilde=start.a2tilde ;
na2=start.na2;
estimafy=start.estimafy; % computes posterior variance
a2tilde_func=start.a2tilde_func ;
%---------------------------------------------

[burnin,thin,waitbar_update]=...
    utils.miscellaneous.parse_arguments(pnames_defaults,varargin{:});

%error checks for burnin and thin
if (burnin<0) || burnin~=round(burnin)
    error('bad burnin value');
end
if (thin<=0)|| thin~=round(thin)
    error('bad thin value');
end

outclass = superiorfloat(start.x0); % single or double

% Put the replicates dimension second.
[npar,nchain] = size(start.x0);

smpl = zeros([npar,nsamples,nchain],outclass);
fsmpl = nan([1,nsamples,nchain],outclass);

% redo this in case the user makes a mistake like giving the negative of
% the log posterior instead of the log posterior directly.
start.f0=nan;
accept_rate = 1;
for idraw = 1-burnin:nsamples*thin
    x1=one_constant_var_posterior_draw();
    
    start.x0 = x1; % preserves x's shape.
    if idraw>0 && mod(idraw,thin)==0; % burnin and thin
        smpl(:,idraw/thin,:) = start.x0;
    end
    if ~isempty(waitbar_update)
        new_messages={
            sprintf('Acceptance rate %0.4f',accept_rate)
            sprintf('thinning factor %0.0f',thin)
            };
        waitbar_update(new_messages{:});
    end
end

    function x1=one_constant_var_posterior_draw()
        
        switch prior_type
            case {'diffuse','jeffrey',1}
                a2tilde.post=estimafy(inv(start.SIGMA),true);
                
                % Posterior of alpha|SIGMA,Data ~ Normal
                alpha2tilde = a2tilde.ols.a + chol(a2tilde.post.V)'*randn(na2,1);% Draw alpha
                
                % Posterior of SIGMA|Data ~ iW(ols.SSE,T-K)
                start.SIGMA = inverse_wishart_draw(a2tilde.ols.SSE,nobs-K);% Draw SIGMA
                
            case {'minnesota',2}
                alpha2tilde = a2tilde.post.a + chol(a2tilde.post.V)'*randn(na2,1); % Draw alpha
                
                % SIGMA in this case is a known matrix, whose form is decided in
                % the prior
                % start.SIGMA=start.SIGMA; % start.SIGMA=a2tilde.prior.SIGMA;
                
            case {'normal_wishart',3}
                % This is the covariance for the posterior density of alpha
                tmp=kron(eye(K),start.SIGMA);
                postValpha2tilde = a2tilde.post.V*a2tilde_func(tmp(start.inv_order,start.inv_order),true);
                
                % Posterior of alpha|SIGMA,Data ~ Normal
                %----------------------------------------
                % in the presence of restrictions, the variance above may
                % not be positive definite and so we replace the cholesky
                % with something else
                [vv,dd]=eig(postValpha2tilde);
                % alpha2tilde = a2tilde.post.a + chol(postValpha2tilde)'*randn(na2,1);
                alpha2tilde = a2tilde.post.a + real(vv*sqrt(dd))*randn(na2,1);
                
                % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
                start.SIGMA = inverse_wishart_draw(a2tilde.post.scale_SIGMA,a2tilde.post.dof_SIGMA);% Draw SIGMA
                
            case {'indep_normal_wishart',4}
                alpha2tilde = a2tilde.post.a + chol(a2tilde.post.V)'*randn(na2,1); % Draw of alpha
                
                ALPHA = start.a2Aprime(alpha2tilde); % Draw of ALPHA
                
                % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
                a2tilde.post.scale_SIGMA = a2tilde.prior.scale_SIGMA + (Y-X*ALPHA.').'*(Y-X*ALPHA.');
                start.SIGMA = inverse_wishart_draw(a2tilde.post.scale_SIGMA,a2tilde.post.dof_SIGMA);% Draw SIGMA
                a2tilde.post=estimafy(inv(start.SIGMA),true);
            otherwise
                error('unknown prior type')
        end
        
        alpha_draw = a_func(alpha2tilde);
        SIGMA_draw = start.SIGMA;
        
        x1=vartools.build_parameter_vector(start.vdata,alpha_draw,SIGMA_draw);
        
    end
end

function A=wishart_draw(S,df)
A = chol(S)'*randn(size(S,1),df);
A = A*A';
end

function A=inverse_wishart_draw(S,df)
A=inv(wishart_draw(inv(S),df));
end

