function [pop,funevals,log_I_init,covar,NumWorkers]=initial_draws(logf,lb,ub,N,penalty,mu,max_attempts)
% initial_draws -- initial draws for populations mcmc algorithms
%
% ::
%
%
%   [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N)
%
%   [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty)
%
%   [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty,mu)
%
% Args:
%
%    - **logf** [function handle]: objective to MINIMIZE
%
%    - **lb** [vector]: lower bound
%
%    - **ub** [vector]: upper bound
%
%    - **N** [integer]: number of parameters to draws
%
%    - **penalty** [numeric|{1e+8}]:
%
%    - **mu** [vector]: start draws
%
%    - **max_attempts** [integer|{50}]: maximum number of attempts for each
%    parameter vector.
%
% Returns:
%    :
%
%    - **pop** [struct]: initial draws, each with fields
%      - **f** [numeric]: value of fitness
%      - **x** [vector]: parameters
%
%    - **funevals** [integer]: number of function evaluations
%
%    - **log_I_init** [0]: initial value of the log marginal data density for
%    a multivariate uniform distribution
%
%    - **covar** [matrix]: variance-covariance matrix of the drawn parameters
%
%    - **NumWorkers** [integer]: number of workers available for parallel
%    processing.
%
% Note:
%
% Example:
%
%    See also:

if nargin<7
    max_attempts=[];
    if nargin<6
        mu=[];
        if nargin<5
            penalty=[];
        end
    end
end
if isempty(max_attempts)
    max_attempts=50;
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

NumWorkers=utils.parallel.get_number_of_workers();

parfor (idraw=1:N,NumWorkers)
    iter=0;
    invalid=true;
    % parfor thinks this is a broadcast variable
    %--------------------------------------------
    local_logf=logf;
    while invalid && iter < max_attempts
        iter=iter+1;
        if iter>1 || any(isnan(theta(:,idraw)))
            theta(:,idraw)=lb+ub_minus_lb.*rand(d,1);
        end
        fx=local_logf(theta(:,idraw));
        funevals(idraw)=funevals(idraw)+1;
        invalid=utils.mcmc.is_violation(fx,penalty);
    end
    if invalid
        error(['no valid parameter after ',int2str(max_attempts),' attempts'])
    end
    pop{idraw}=struct('f',fx,'x',theta(:,idraw));
end
pop=cell2mat(pop);
pop=utils.mcmc.sort_population(pop);
x=[pop.x];
covar=cov(x.');
end
