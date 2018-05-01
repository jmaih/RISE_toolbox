function [x0,retcode,tau]=transpose_free_quasi_minimum_residual(...
A,... % coefficient matrix
b,... % right hand side
x0,... % initial guess
tol,... % tolerance level
MaxIter,... % maximum number of iterations
verbose) % flag for printing progress or not

% transpose_free_quasi_minimum_residual attempts to solve Ax=b
%
% ::
%
%   [x0,retcode,tau]=transpose_free_quasi_minimum_residual(A,b,x0)
%   [x0,retcode,tau]=transpose_free_quasi_minimum_residual(A,b,x0,tol)
%   [x0,retcode,tau]=transpose_free_quasi_minimum_residual(A,b,x0,tol,MaxIter)
%   [x0,retcode,tau]=transpose_free_quasi_minimum_residual(A,b,x0,tol,MaxIter,verbose)
%
% Args:
%    - A :
%    - b :
%    - x0:
%    - tol:
%    - MaxIter:
%    - verbose :
%
% Returns:
%    :
%    - x0 :
%    - retcode :
%    - tau :
%
% Note:
%
% Example:
%
%    See also:


% attempts to solve Ax=b

%    retcode: 0 TFQMR converged to the desired tolerance TOL within MAXIT iterations.
%    retcode: 201 TFQMR iterated MAXIT times without converging. This may
%             still be a good solution. 
%    retcode: 202 TFQMR delivered Nans: probably no solution.

% examples:
% m=7;
% AA=rand(m);
% bb=rand(m,1);
% x0=rand(m,1);
% % first input could be a function
% xopt=transpose_free_quasi_minimum_residual(@(q)AA*q,bb,x0,[]);
% first input could be a square matrix
% xopt2=transpose_free_quasi_minimum_residual(AA,bb,x0);
% xopt3=transpose_free_quasi_minimum_residual(AA,bb,x0,sqrt(eps));
% xopt4=transpose_free_quasi_minimum_residual(AA,bb,x0,sqrt(eps),20);
% xopt5=transpose_free_quasi_minimum_residual(AA,bb,x0,sqrt(eps),20,true));
% norm(AA*[xopt,xopt2,xopt3,xopt4,xopt5]-bb(:,ones(1,5)))

%   References:
%      1. R. Freund, A transpose-free quasi-minimal residual algorithm for
%         non-Hermitian linear systems, SIAM J. Sci. Comp., 14 (1993),
%         pp. 470--482.
%      2. C. T. Kelley, Iterative Methods for Linear and Nonlinear equations, 2nd ed.
%         SIAM, 1995, Philadelphia. pp. 55
%      2. Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd ed.
%         SIAM, 2003, Philadelphia.


switch class(A)
    case {'inline','function_handle'}
        matrix_times_vector=A;
    case {'double','float'}
        matrix_times_vector=@(x)A*x;
    case 'char'
        matrix_times_vector=str2func(A);
    otherwise
        error([mfilename,':: A should be a matrix, a function handle, an inline function or a string representing a function'])
end
n=size(b,1);

assert(nargin>=3 && nargin<=6,'number of arguments should be between 3 and 6')
if nargin<6
    verbose=false;
    if nargin<5
        MaxIter=n;
        if nargin<4
            tol=sqrt(eps);
        end
    end
end
n2b=norm(b);
e_norm_b=tol*n2b;

retcode=0;
if isempty(x0)
    x0=zeros(n,1);
end
if e_norm_b==0
    x0=zeros(n,1);
    tau=0;
else
    stagnated=0;
    max_stagnations=3;
    r0=b;
    if any(x0)
        r0=r0-matrix_times_vector(x0);
    end
    w=r0; % <--- w=r0(:,[1,1]);
    y=w;
    Ay=matrix_times_vector(y); % <--- v=matrix_times_vector(y(:,1));
    v=Ay;
    tau=norm(r0);
    rho=r0'*r0;
    theta=0;
    eta=0;
    iter=0;
    d=0;
    [converged,finished]=assess_convergence(0);
    while ~finished
        iter=iter+1;
        sigma=r0'*v;
        alpha=rho/sigma;
        for jj=1:2
            if ~converged
                m=jj;
                if jj==2
                    y=y-alpha*v; 
                    Ay=matrix_times_vector(y); 
                end
                d=y+(theta^2*eta/alpha)*d; 
                w=w-alpha*Ay; 
                theta=norm(w)/tau; 
                c=1/sqrt(1+theta^2);
                tau=tau*theta*c; % == norm(r)/norm(b)
                eta=c^2*alpha;
                xtmp=x0+eta*d;
                stag_flag=norm(xtmp-x0)<1e-8;
                stagnated=stag_flag*(stagnated+stag_flag);
                x0=xtmp;
                [converged,finished]=assess_convergence(2*(iter-1)+m);
            end
        end
        if verbose
            fprintf(1,'%s %6.0d %s %6.4d\n','iter',iter,'tau',tau);
        end
        if ~finished
            rho_new=r0'*w; 
            beta=rho_new/rho;
            % partial update of v
            v=beta*(Ay+beta*v);
            y=w+beta*y; 
            Ay=matrix_times_vector(y); 
            v=Ay+v; 
            rho=rho_new;
        end
    end
    if isnan(tau)
        retcode=202;
    elseif iter>=MaxIter
        if (tau*sqrt(m+1)>e_norm_b)
            retcode=201;
        end
    elseif stagnated>=max_stagnations
        %     disp([mfilename,':: stagnated before converging properly. May be useful to use a higher tolerance level'])
    end
end

    function [conv,finished]=assess_convergence(m)
        conv= tau*sqrt(m+1)<=e_norm_b;
        finished=conv||...
            iter>=MaxIter||...
            stagnated>=max_stagnations;
    end
end