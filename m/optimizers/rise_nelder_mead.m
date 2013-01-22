function [xbest,fbest,exitflag,output] = rise_nelder_mead(objective,xbest,lb,ub,options,varargin)
% simplex algorithm with nonlinear constraints

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

default_options = struct('MaxIter',100000,...
    'MaxFunEvals',150000,'TolX',1e-4,'TolFun',1e-4,'nonlcon',@(z)0);

% If just 'defaults' passed in, return the default options in X
if nargin==0
    xbest=default_options;
    return
end

if nargin<5
    options=[];
end

ff=fieldnames(default_options);
for ifield=1:numel(ff)
    v=ff{ifield};
    if isfield(options,v)
        default_options.(v)=options.(v);
    end
end

MaxFunEvals=default_options.MaxFunEvals;% : iteration limit in terms of function evaluations
nonlcon=default_options.nonlcon;
MaxIter=default_options.MaxIter;
TolX=default_options.TolX;
TolFun=default_options.TolFun;

%         [objective,xbest,options] = separateOptimStruct(objective);
n = numel(xbest);
funcCount = 0;
iterations = 0;

header = ' Iteration   Func-count     min f(xbest)         Procedure';

% Initialize parameters
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;

% Set up a simplex near the initial guess.

pop = wrapper();
pop(1) = wrapper(xbest);
how = '';
% Initial simplex setup continues later

disp(' ')
disp(header)
fprintf(' %5.0f        %5.0f     %12.6g         %s\n', iterations, funcCount, pop(1).f, how);

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;             % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;      % Even smaller delta for zero elements of xbest
for j = 1:n
    y = xbest;
    if y(j) ~= 0
        y(j) = (1 + usual_delta)*y(j);
    else
        y(j) = zero_term_delta;
    end
    pop(1,j+1) = wrapper(y);
end

% sort so pop(1,:) has the lowest function value
pop = sort_population(pop);

how = 'initial simplex';
iterations = iterations + 1;
fprintf(' %5.0f        %5.0f     %12.6g         %s\n', iterations, funcCount, pop(1).f, how)

% Main algorithm: iterate until
% (a) the maximum coordinate difference between the current best point and the
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v2-v1||,...,||pop(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations
% are exceeded
while funcCount < MaxFunEvals && iterations < MaxIter
    if max(abs(pop(1).f-[pop(2:n+1).f])) <= max(TolFun,10*eps(pop(1).f)) && ...
            max(max(abs(bsxfun(@minus,[pop(2:n+1).x],pop(1).x)))) <= max(TolX,10*eps(max(pop(1).x)))
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the n (NOT n+1) best points
    xbar = sum([pop(:,1:n).x],2)/n;
    xfr=newpoint('reflect');
    
    choice=compare_individuals(xfr,pop(1));
    if choice==1
        % Calculate the expansion point
        xfe=newpoint('expand');
        choice=compare_individuals(xfe,pop(1));
        if choice==1
            pop(end) = xfe;
            how = 'expand';
        else
            pop(end) = xfr;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        choice=compare_individuals(xfr,pop(n));
        if choice==1
            pop(end) = xfr;
            how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            choice=compare_individuals(xfr,pop(end));
            if choice==1
                % Perform an outside contraction
                xfc=newpoint('contract outside');
                choice=compare_individuals(xfc,xfr);
                if choice==1
                    pop(end) = xfc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xfcc=newpoint('contract inside');
                 choice=compare_individuals(xfcc,pop(end));
                if choice==1
                    pop(end) = xfcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                pop=newpoint('shrink');
            end
        end
    end
    pop = sort_population(pop);
    iterations = iterations + 1;
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', iterations, funcCount, pop(1).f, how)
end   % while

xbest = pop(1).x;
fbest =  pop(1).f;

output.iterations = iterations;
output.funcCount = funcCount;
output.algorithm = mfilename;

if funcCount >= MaxFunEvals || iterations >= MaxIter
    exitflag = 0;
else
    exitflag = 1;
end

    function indiv=newpoint(type)
        switch type
            case 'shrink'
                for ipoint=2:n+1
                    xpoint=pop(1).x+sigma*(pop(ipoint).x-pop(1).x);
                    pop(ipoint)= wrapper(xpoint);
                end
                indiv= pop;
                return
            case 'contract outside'
                coef=psi*rho;
            case 'contract inside'
                coef=-psi;
            case 'reflect'
                coef=rho;
            case 'expand'
                coef=rho*chi;
        end
        xbest = xbar+ coef*(xbar-pop(end).x);
        indiv = wrapper(xbest);
    end

    function this=wrapper(bird)
        if nargin==0
            this=evaluate_individual();
        else
            this=evaluate_individual(bird,objective,lb,ub,nonlcon,varargin{:});
            funcCount=funcCount+1;
        end
    end

end

