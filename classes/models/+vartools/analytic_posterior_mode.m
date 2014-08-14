function [a_post,iV_post,resid]=analytic_posterior_mode(X,y,a_prior,Va_prior)

%     a_prior=a2tilde_func(a_prior);
%     Va_prior=a2tilde_func(va_prior,true);
nparam=numel(a_prior);
iVa_prior=Va_prior\eye(nparam);

% ols
%----
a_ols=X\y;
iV_ols=(X'*X);

% posterior
%----------
iV_post=(iV_ols+iVa_prior);
a_post=iV_post\(iV_ols*a_ols+iVa_prior*a_prior);

% residuals
resid=y-X*a_post;
end