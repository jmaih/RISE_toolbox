function [x1]=one_constant_var_posterior_draw(X,Y,...
    post_a2tilde,ols_a2tilde,prior_a2tilde,...
    a_func,K,endo_nbr,nobs,prior_type)

persistent SIGMA a2A

if isempty(SIGMA)
    SIGMA = ols_a2tilde.SIGMA; % This is the single draw from the posterior of SIGMA
    % NB: We transpose here so that we can do the computations as in the basic
    % regression analysis !!!!
    %-------------------------------------------------------------------------
    a2A=@(x)transpose(reshape(a_func(x),endo_nbr,K));
end
XpX=X'*X;

A_OLS=a2A(ols_a2tilde.a);

if any(strcmp(prior_type,{'normal_wishart','indep_normal_wishart'}))
    % Hyperparameters on inv(SIGMA) ~ W(prior.dof_SIGMA,inv(prior.scale_SIGMA))
    prior_a2tilde.dof_SIGMA = endo_nbr+1;         %<---- prior Degrees of Freedom (DoF) of SIGMA
    prior_a2tilde.scale_SIGMA = eye(endo_nbr);    %<---- prior scale of SIGMA
    % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
    post_a2tilde.dof_SIGMA = nobs + prior_a2tilde.dof_SIGMA;
    if strcmp(prior_type,'normal_wishart')
        % we invert and then apply the function:probably the simplest thing
        % to do
        %------------------------------------------------------------------
        iVa=a_func(inv(prior_a2tilde.V),true); 
        post_a2tilde.scale_SIGMA = ols_a2tilde.SSE + prior_a2tilde.scale_SIGMA + ...
            A_OLS'*XpX*A_OLS + ...
            A_prior'*iVa*A_prior - ...
            A_post'*(iVa + XpX)*A_post;
    end
end

switch prior_type
    case 'diffuse' % prior == 1
        % Posterior of alpha|SIGMA,Data ~ Normal
        alpha2tilde = ols_a2tilde.a + chol(post_a2tilde.V)'*randn(na2,1);% Draw alpha
        
        % Posterior of SIGMA|Data ~ iW(ols.SSE,T-K)
        SIGMA = inverse_wishart_draw(ols_a2tilde.SSE,nobs-K);% Draw SIGMA
        
    case 'minnesota' % prior == 2
        alpha2tilde = post_a2tilde.a + chol(post_a2tilde.V)'*randn(na2,1); % Draw alpha
        
        % SIGMA in this case is a known matrix, whose form is decided in
        % the prior
        SIGMA=prior_a2tilde.SIGMA;
        
    case 'normal_wishart' % prior == 3
        % This is the covariance for the posterior density of alpha
        COV = kron(SIGMA,post_a2tilde.V);
        
        % Posterior of alpha|SIGMA,Data ~ Normal
        alpha2tilde = post_a2tilde.a + chol(COV)'*randn(na2,1);  % Draw alpha
        
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
        SIGMA = inverse_wishart_draw(post_a2tilde.scale_SIGMA,post_a2tilde.dof_SIGMA);% Draw SIGMA
        
    case 'indep_normal_wishart' % prior == 4
        alpha2tilde = post_a2tilde.a + chol(post_a2tilde.V)'*randn(na2,1); % Draw of alpha
        
        ALPHA = a2A(alpha2tilde); % Draw of ALPHA
        
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(post.scale_SIGMA),post.dof_SIGMA)
        post_a2tilde.scale_SIGMA = prior_a2tilde.scale_SIGMA + (Y-X*ALPHA)'*(Y-X*ALPHA);
        SIGMA = inverse_wishart_draw(post_a2tilde.scale_SIGMA,post_a2tilde.dof_SIGMA);% Draw SIGMA
    otherwise
        error('unknown prior type')
end

alpha_draw = a_func(alpha2tilde);
SIGMA_draw = SIGMA;

x1=vartools.build_parameter_vector(vdata,alpha_draw,SIGMA_draw);

end

function A=wishart_draw(S,df)
A = chol(S)'*randn(size(S,1),df);
A = A*A';
end

function A=inverse_wishart_draw(S,df)
A=inv(wishart_draw(inv(S),df));
end
