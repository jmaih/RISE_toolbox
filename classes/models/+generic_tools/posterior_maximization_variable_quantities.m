function post_max=posterior_maximization_variable_quantities(post_max,H,flag)
% posterior_maximization_variable_quantities -- recomputes the quantities
% that depend of the Hessian
%
% Syntax
% -------
% ::
%
%   post_max=posterior_maximization_variable_quantities(post_max,H)
%
%   post_max=posterior_maximization_variable_quantities(post_max,H,flag)
%
% Inputs
% -------
%
% - **post_max** [struct]: more specifically the content of
% obj.estimation.posterior_maximization
%
% - **H** [d x d x k array]: Array of hessians. Each page represents a
% hessian.
%
% - **flag** [{true}|false|numeric]: indicator for the type of hessian to
% use. If true, the hessian returned by the optimizer is used. If false,
% then the hessian computed numerically is used.
%
% Outputs
% --------
%
% - **post_max** [struct]: containing
%   - **SD** [vector]: standard deviations of the parameters
%   - **log_marginal_data_density_laplace** [numeric]: laplace
%   approximation of the log marginal data density
%   - **vcov** [matrix]: variance-covariance of estimated parameters
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% flag = true --> optimizer hessian
% flag = false --> numerical hessian
if nargin<2
    flag=true;
end
if islogical(flag)
    flag=double(flag);
end
if flag==0
    flag=2;
end

pages=1:size(H,3);
if isempty(pages)
    error('no valid hessian entered')
end
if pages(end)>1
    if ~any(flag-pages==0)
        error('index for the choice of hessian exceeds third dimension of H')
    end
end
if is_valid_hessian(H(:,:,flag))
    H=H(:,:,flag);
else
    warning('optimizer hessian unavailable or invalid, reverting to numerical hessian')
    H=H(:,:,2);
end

d=size(H,1);
Hinv=H\eye(d);
SD=sqrt(diag(Hinv));
post_max.mode_stdev=SD; %
post_max.log_marginal_data_density_laplace=...
    utils.marginal_data_density.laplace(post_max.log_post,Hinv);
post_max.vcov=Hinv; %

    function flag=is_valid_hessian(H)
        flag=~any(isnan(vec(H)));
    end

end