function [Results,accept_rate,start] = constant_bvar_sampler(start,options)
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

if nargin<2
    options=[];
end
if isempty(options)
    options=struct();
end

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x) && isreal(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
defaults={ % arg_names -- defaults -- checks -- error_msg
    'burnin',0,@(x)num_fin_int(x),'burnin should be an integer in [0,inf)'
    'N',20000,@(x)num_fin_int(x) && x>0 ,'N should be a strictly positive integer'
    'thin',1,@(x)num_fin_int(x) && x>=1,'thin must be >=1'
    'nchain',1,@(x)num_fin_int(x) && x>0,'nchain(# chains) should be a strictly positive integer'
    'MaxTime',inf,@(x)num_fin(x) && x>0,'MaxTime must be a positive scalar'
    };
if nargin==0
    Results=cell2struct(defaults(:,2),defaults(:,1),1);
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

options=utils.miscellaneous.parse_arguments(defaults,options);

MaxTime=options.MaxTime;
nchain=options.nchain;
burnin=options.burnin;
thin=options.thin;
N=options.N;

% Put the replicates dimension second.
[npar,nchain] = size(start.x0);

smpl = zeros(npar,N,nchain);

% redo this in case the user makes a mistake like giving the negative of
% the log posterior instead of the log posterior directly.
start.f0=nan;

obj=struct();
obj.funcCount=sum(funevals);
obj.iterations=0;
obj.start_time=clock;
obj.MaxTime=options.MaxTime;
obj.MaxFunEvals=inf;

accept_rate = 1;
idraw=-burnin;
total_draws=N*thin+burnin;
obj.MaxIter=total_draws;
utils.optim.manual_stopping();
stopflag=utils.optim.check_convergence(obj);
wtbh=waitbar(0,'please wait...','Name','MH sampling');
accept_ratio=[];
while isempty(stopflag)
    idraw=idraw+1;
    obj.iterations=idraw+burnin;
    x1=one_constant_var_posterior_draw();
    
    start.x0 = x1; % preserves x's shape.
    if idraw>0 && mod(idraw,thin)==0; % burnin and thin
        smpl(:,idraw/thin,:) = start.x0;
    end
    % display progress
    %------------------
    if mod(obj.iterations,100)==0
        fprintf(1,'iter # %0.0f, Acceptance rate %0.4f, thinning factor %0.0f\n',...
            obj.iterations,accept_rate,thin);
    end
    stopflag=utils.optim.check_convergence(obj);
    waitbar_updating()
end
delete(wtbh)

% pre-allocate
%--------------
Results=struct();
Results.pop=struct('f',nan,'x',nan);
Results.pop=Results.pop(ones(1,nchain),ones(1,N));
for isample=1:N
    for ichain=1:nchain
        Results.pop(ichain,isample).x=smpl(:,isample,ichain);
    end
end
Results.m=mean(smpl,2);

    function waitbar_updating()
        x=obj.iterations/total_draws;
        waitbar(x,wtbh,...
            {
            sprintf('bestf %s',num2str(bestf))
            sprintf('acceptance rate %s',num2str(100*accept_ratio))
            }...
            )
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

