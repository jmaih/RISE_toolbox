function post_max=posterior_maximization_variable_quantities(post_max,...
a_func)
% posterior_maximization_variable_quantities -- recomputes the quantities
% that depend on the Hessian
%
% ::
%
%
%   post_max=posterior_maximization_variable_quantities(post_max,H)
%
%   post_max=posterior_maximization_variable_quantities(post_max,H,flag)
%
% Args:
%
%    - **post_max** [struct]: more specifically the content of
%    obj.estimation.posterior_maximization
%
%    - **a_func** [function_handle]: function that inflates x and Vx under
%    linear restrictions.
%
% Returns:
%    :
%
%    - **post_max** [struct]: containing
%      - **SD** [vector]: standard deviations of the parameters
%      - **log_marginal_data_density_laplace** [numeric]: laplace
%      approximation of the log marginal data density
%      - **vcov** [matrix]: variance-covariance of estimated parameters
%
% Note:
%
% Example:
%
%    See also:

H=post_max.hessian;

d=size(H,1);

flag=is_valid_hessian(H);

Hinv=nan(size(H));

if flag
    
    Hinv=H\eye(d);
    
else
    
    disp('-----------------------------------------------------')
    disp([mfilename,upper(':: gentle warning: the hessian is invalid')])
    disp('-----------------------------------------------------')
    
end

% compute marginal data density based on the short hessian
%----------------------------------------------------------

post_max.log_marginal_data_density_laplace=...
    utils.marginal_data_density.laplace_mdd(post_max.log_post,Hinv);

% inflate the covariance under linear restrictions
%--------------------------------------------------
try
    
    post_max.vcov=a_func(Hinv,true); %
    
catch
    % under constant-parameter RFVAR we do not compute the hessian at
    % posterior maximization and if there are parameter restrictions, the
    % step above will not work
    post_max.vcov=Hinv; %
    
end

% compute the standard deviations based on the inflated covariance
%------------------------------------------------------------------
SD=sqrt(diag(post_max.vcov));

post_max.mode_stdev=SD; %

    function flag=is_valid_hessian(H)
        
        flag=~any(isnan(vec(H)));
        
    end

end