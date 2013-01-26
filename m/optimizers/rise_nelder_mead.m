function [xbest,fbest,exitflag,output] = rise_nelder_mead(objective,x0,lb,ub,options,varargin)
% simplex algorithm with nonlinear constraints

%   References:
%   Lagarias, J. C., J. A. Reeds, M. H. Wright,  P. E. Wright(1998):
%   "Convergence Properties of the Nelder-Mead Simplex Method in Low
%    Dimensions", SIAM Journal of Optimization, 9(1),p.112-147.
%   Gao, F., L. Lixing(2010): "Implementing the Nelder-mead simplex
%   algorithm with adaptive parameters", Comput Optim Appl.

default_options = struct('MaxIter',100000,...
    'MaxFunEvals',150000,'TolX',1e-4,'TolFun',1e-4,...
    'nonlcon',@(z)0,...
    'simplex_method',[]);

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
simplex_method=default_options.simplex_method;
if isempty(simplex_method)
    simplex_method='';
end

%         [objective,xbest,options] = separateOptimStruct(objective);
npar = numel(x0);
funcCount = 0;
iterations = 0;

header = ' Iteration   Func-count     min f(x)         Procedure';

% Initialize parameters
% reflect_param = 1; % rho
% expand_param = 2; % chi
% shrink_param = 0.5; % psi
% sigma = 0.5;

reflect_param = 1; % rho
expand_param = 1+2/npar; % chi
shrink_param = max(0.75-1/(2*npar),.5); % psi
sigma = max(1-1/npar,0.5);

% Set up a simplex near the initial guess.

Simplex = wrapper();
construct_simplex();
how = '';
% Initial simplex setup continues later

disp(' ')
disp(header)
fprintf(' %5.0f        %5.0f     %12.6g         %s\npar', iterations, funcCount, Simplex(1).f, how);

% sort so Simplex(1,:) has the lowest function value
Simplex = sort_population(Simplex);

how = 'initial simplex';
iterations = iterations + 1;
fprintf(' %5.0f        %5.0f     %12.6g         %s\npar', iterations, funcCount, Simplex(1).f, how)

% Main algorithm: iterate until
% (a) the maximum coordinate difference between the current best point and the
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v2-v1||,...,||Simplex(npar+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations
% are exceeded
while funcCount < MaxFunEvals && iterations < MaxIter
    if max(abs(Simplex(1).f-[Simplex(2:npar+1).f])) <= max(TolFun,10*eps(Simplex(1).f)) && ...
            max(max(abs(bsxfun(@minus,[Simplex(2:npar+1).x],Simplex(1).x)))) <= max(TolX,10*eps(max(Simplex(1).x)))
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the npar (NOT npar+1) best points
    xbar = sum([Simplex(:,1:npar).x],2)/npar;
    xf_r=newpoint('reflect');
    
    choice=compare_individuals(xf_r,Simplex(1));
    if choice==1
        % Calculate the expansion point
        xf_e=newpoint('expand');
        choice=compare_individuals(xf_e,Simplex(1));
        if choice==1
            Simplex(end) = xf_e;
            how = 'expand';
        else
            Simplex(end) = xf_r;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        choice=compare_individuals(xf_r,Simplex(npar));
        if choice==1
            Simplex(end) = xf_r;
            how = 'reflect';
        else % fxr >= fv(:,npar)
            % Perform contraction
            choice=compare_individuals(xf_r,Simplex(end));
            if choice==1
                % Perform an outside contraction
                xf_oc=newpoint('contract outside');
                choice=compare_individuals(xf_oc,xf_r);
                if choice==1
                    Simplex(end) = xf_oc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xf_ic=newpoint('contract inside');
                choice=compare_individuals(xf_ic,Simplex(end));
                if choice==1
                    Simplex(end) = xf_ic;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                Simplex=newpoint('shrink');
            end
        end
    end
    Simplex = sort_population(Simplex);
    iterations = iterations + 1;
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', iterations, funcCount, Simplex(1).f, how)
end   % while

xbest = Simplex(1).x;
fbest =  Simplex(1).f;

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
            case 'contract outside'
                coef=shrink_param*reflect_param;
            case 'contract inside'
                coef=-shrink_param;
            case 'reflect'
                coef=reflect_param;
            case 'expand'
                coef=reflect_param*expand_param;
            case 'shrink'
                for ipoint=2:npar+1
                    xpoint=Simplex(1).x+sigma*(Simplex(ipoint).x-Simplex(1).x);
                    Simplex(ipoint)= wrapper(xpoint);
                end
                indiv= Simplex;
        end
        if ~strcmp(type,'shrink')
            xbest = xbar+ coef*(xbar-Simplex(end).x);
            indiv = wrapper(xbest);
        end
    end

    function this=wrapper(bird)
        if nargin==0
            this=evaluate_individual();
        else
            this=evaluate_individual(bird,objective,lb,ub,nonlcon,varargin{:});
            funcCount=funcCount+1;
        end
    end

    function construct_simplex()
        % Continue setting up the initial simplex.
        % Following improvement suggested by L.Pfeffer at Stanford
        usual_delta = 0.05;             % 5 percent deltas for non-zero terms
        zero_term_delta = 0.00025;      % Even smaller delta for zero elements of xbest
        Simplex(1) = wrapper(x0);
        for j = 1:npar
            switch simplex_method
                case ''
                    y = x0;
                    if y(j) ~= 0
                        y(j) = (1 + usual_delta)*y(j);
                    else
                        y(j) = zero_term_delta;
                    end
                case 'random'
                    y=lb+(ub-lb).*rand(npar,1);
                otherwise
                    error(['unknown simplex construction method ',simplex_method])
            end
            Simplex(1,j+1) = wrapper(y);
        end
    end

end

