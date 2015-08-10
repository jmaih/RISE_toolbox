function [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty,mu)
% initial_draws -- initial draws for populations mcmc algorithms
%
% Syntax
% -------
% ::
%
%   [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N)
%
%   [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty)
%
%   [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty,mu)
%
% Inputs
% -------
%
% - **logf** [function handle]: objective to MINIMIZE  
%
% - **lb** [vector]: lower bound  
%
% - **ub** [vector]: upper bound  
%
% - **N** [integer]: number of parameters to draws  
%
% - **penalty** [numeric|{1e+8}]:  
%
% - **mu** [vector]: start draws  
%
% Outputs
% --------
%
% - **pop** [struct]: initial draws, each with fields  
%   - **f** [numeric]: value of fitness
%   - **x** [vector]: parameters  
%
% - **funevals** [integer]: number of function evaluations  
%
% - **log_I_init** [0]: initial value of the log marginal data density for
% a multivariate uniform distribution
%
% - **covar** [matrix]: variance-covariance matrix of the drawn parameters  
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if nargin<6
    mu=[];
    if nargin<5
        penalty=[];
    end
end
if isempty(penalty)
    penalty=1e+8;
end
% this has to be a probability density and not a Kernel and so it
% sums to 1. Do initialization for each core and for the general...
I_init=1;
log_I_init=log(I_init);

d=size(lb,1);
n0=size(mu,2);

ub_minus_lb=ub-lb;
% initial vector is always part of the game
%------------------------------------------
theta=[mu,nan(d,N-n0)];
pop=cell(1,N);
funevals=zeros(1,N);
for idraw=1:N
    iter=0;
    invalid=true;
    while invalid && iter < maxiter
        iter=iter+1;
        if iter>1 || ~isnan(theta(:,idraw))
            theta(:,idraw)=lb+ub_minus_lb.*rand(d,1);
        end
        fx=logf(theta(:,idraw));
        funevals=funevals+1;
        invalid=utils.mcmc.is_violation(fx,penalty);
    end
    if invalid
        error(['no valid parameter after ',int2str(maxiter),' attempts'])
    end
    pop{idraw}=struct('f',fx,'x',theta(:,idraw));
end
pop=cell2mat(pop);
pop=utils.mcmc.sort_population(pop);
x=[pop.x];
covar=cov(x.');
end
