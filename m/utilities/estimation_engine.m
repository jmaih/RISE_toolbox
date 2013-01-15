function [xbest,fbest,H,issue]=estimation_engine(PROBLEM,hessian_type)

opt.optimset=optimset('Display','iter',...%[ off | iter | iter-detailed | notify | notify-detailed | final | final-detailed ]
    'MaxFunEvals',inf,...% [ positive scalar ]
    'MaxIter',1000,...%: [ positive scalar ]
    'TolFun',sqrt(eps),...%: [ positive scalar ]
    'MaxNodes',20,...%: [ positive scalar | {1000*numberOfVariables} ]
    'UseParallel','always',...%: [ always | {never} ]
    'MaxTime',3600);%: [ positive scalar | {7200} ]
opt.optimizer='fmincon'; % default is fmincon
opt.hessian_type='fd';
% 'Algorithm','interior-point',...% [ active-set | interior-point | levenberg-marquardt | sqp | ...
% I observed the following problem with the interior point algorithm (IP): I estimated a model using my abc
% routines and then wanted to sharpen the results using fmincon. fmincon with the IP
% could not even replicate the function at the starting point (my guess is that this is
% because some parameters were lying against their boundaries). More exactly, the function was given the
% the value of the penalty. When I turned off the interior point, everything was fine. Maybe all this is
% consistent with the definition of interior point.
% further testing with wrapper function fminsearch_bnd also gave the wrong starting value function. This
% points in the direction of a penalty arising with the parameters not being in the interior.
% I probably should revisit the way I impose constraints in abc

if nargin==0
    if nargout>1
        error([mfilename,':: number of output arguments cannot exceed 1 when there are no input arguments'])
    end
    xbest=opt;
    return
end
OtherProblemFields={'Aineq','bineq','Aeq','beq','nonlcon'};
if isempty(PROBLEM.solver)
    PROBLEM.solver=opt.optimizer;
end
for ii=1:numel(OtherProblemFields)
    if ~isfield(PROBLEM,OtherProblemFields{ii});
        PROBLEM.(OtherProblemFields{ii})=[];
    end
end
if isnumeric(PROBLEM.solver)
    switch PROBLEM.solver
        case 0
        case 1
            PROBLEM.solver='fmincon';
        case 2
            PROBLEM.solver='fminunc';
        case 3
        case 4
        case 5
        case 6
        case 7
            PROBLEM.solver='fminsearch';
        otherwise
    end
end
optimizer=PROBLEM.solver;
fh=PROBLEM.objective;
x0=PROBLEM.x0;
lb=PROBLEM.lb;
ub=PROBLEM.ub;
optim_opt=PROBLEM.options;

if isempty(hessian_type)
    hessian_type=opt.hessian_type;
end
options=opt.optimset;
options=mysetfield(options,optim_opt);

if isempty(hessian_type)
    hessian_type='fd';
end

Nsim=size(x0,2);
xfinal=cell(1,Nsim);
ffinal=cell(1,Nsim);
exitflag=cell(1,Nsim);

fh_handle=@fh_wrapper;
problem_maker=@make_problem;
if iscell(optimizer)
    vargs=optimizer(2:end);
    optimizer=optimizer{1};
    for ii=1:Nsim
        [xfinal{ii},ffinal{ii},exitflag{ii}]=optimizer(fh_handle,x0(:,ii),lb,ub,options,vargs{:});
    end
elseif isa(optimizer,'function_handle')
    for ii=1:Nsim
        [xfinal{ii},ffinal{ii},exitflag{ii}]=optimizer(problem_maker(ii));
    end
else
    switch optimizer
        case 0 % just compute log-likelihood/posterior
            for ii=1:Nsim
                ffinal{ii}=fh(x0(:,ii));
                xfinal{ii}=x0(:,ii);
                exitflag{ii}=1;
            end
        case {1,'fmincon'}
            for ii=1:Nsim
                [xfinal{ii},ffinal{ii},exitflag{ii}] = fmincon(problem_maker(ii));
            end
        case {2,'fminunc'}
            for ii=1:Nsim
                [xfinal{ii},ffinal{ii},exitflag{ii}]=fminunc_bnd(fh_handle,...
                    x0(:,ii),lb,ub,options);
            end
        case {4,5,6}
            error([mfilename,':: cases 4, 5 and 6 are not yet implemented'])
        case {7,'fminsearch'} %
            for ii=1:Nsim
                [xfinal{ii},ffinal{ii},exitflag{ii}]=fminsearch_bnd(fh_handle,...
                    x0(:,ii),lb,ub,options);
            end
        otherwise
            if ischar(optimizer)
                optimizer=str2func(optimizer);
            end
            for ii=1:Nsim
                [xfinal{ii},ffinal{ii},exitflag{ii}]=optimizer(problem_maker(ii));
            end
    end
end

% select the overall best
[~,order]=sort(cell2mat(ffinal));
best=order(1);
xbest=xfinal{best};
fbest=ffinal{best};
TranslateOptimization(exitflag{best});

% compute Hessian
issue='';
switch lower(hessian_type)
    case 'fd'
        H = finite_difference_hessian(fh_handle,xbest);
    case 'opg'
        H = outer_product_hessian(fh_handle,xbest);
        if any(any(isnan(H)))||any(any(isinf(H)))
            issue='OPG unstable and inaccurate for calculation of Hessian, switched to finite differences';
            warning([mfilename,':: ',issue]) %#ok<WNTAG>
            warning([mfilename,':: OPG unstable for calculation of Hessian, switching to finite differences']) %#ok<WNTAG>
            H = finite_difference_hessian(fh_handle,xbest);
        end
    otherwise
        issue=['unknow hessian option ',hessian_type,' using finite differences'];
        warning([mfilename,':: ',issue]) %#ok<WNTAG>
        H = finite_difference_hessian(fh_handle,xbest);
end

% % %     function [constr,constreq]=my_nonlinear_restriction(x)
% % %         constr=nonlcon(x);
% % %         constreq=[];
% % %     end

    function f=fh_wrapper(x)
        f=fh(x);
    end

    function pb=make_problem(index)
        pb=PROBLEM;
        pb.x0=x0(:,index);
    end
end

function TranslateOptimization(exitflag)
switch exitflag
    case 1
        disp('First order optimality conditions satisfied.')
    case 0
        disp('Too many function evaluations or iterations.')
    case -1
        disp('Stopped by output/plot function.')
    case -2
        disp('No feasible point found.')
    case 2
        %     Trust-region-reflective and interior-point:
        disp('Change in X too small.')
    case 3
        %     Trust-region-reflective:
        disp('Change in objective function too small.')
    case 4
        %     Active-set only:
        disp('Computed search direction too small.')
    case 5
        disp('Predicted change in objective function too small.')
    case -3
        %     Interior-point:
        disp('Problem seems unbounded.')
    otherwise
        % do nothing
end
end

function [x,fval,exitflag]=fminsearch_bnd(fun,x0,lb,ub,options)
bc=bound_class([lb,ub]);
xt=transform_parameters(x0,bc,[lb,ub]);
[xt,fval,exitflag]=fminsearch(@wrapper,xt,options,fun,bc,lb,ub);
x=untransform_parameters(xt,bc,[lb,ub]);
end

function [x,fval,exitflag]=fminunc_bnd(fun,x0,lb,ub,options)
bc=bound_class([lb,ub]);
xt=transform_parameters(x0,bc,[lb,ub]);
[xt,fval,exitflag]=fminunc(@wrapper,xt,options,fun,bc,lb,ub);
x=untransform_parameters(xt,bc,[lb,ub]);
end

function ff=wrapper(xt,fun,bc,lb,ub)
xu=untransform_parameters(xt,bc,[lb,ub]);
ff=fun(xu);
end

function bc=bound_class(bounds)
bc=isinf(bounds(:,1))+2*isinf(bounds(:,2));
end

function xt=transform_parameters(xu,bc,bounds)
xt=xu;
% unbounded, id=bc==3; xt(id)=xu(id);
% do nothing
% lower bound and upper bound
id=bc==0;
xt(id)=log((xu(id)-bounds(id,1))./(bounds(id,2)-xu(id)));
% lower bound only
id=bc==2;
xt(id)=log(xu(id)-bounds(id,1));
% upper bound only
id=bc==1;
xt(id)=log(bounds(id,2)-xu(id));
end

function xu=untransform_parameters(xt,bc,bounds)
xu=xt;
% unbounded, id=bc==3; xt(id)=xu(id);
% do nothing
% lower bound and upper bound
id=bc==0;
xu(id)=(bounds(id,1)+bounds(id,2).*exp(xt(id)))./(1+exp(xt(id)));
% lower bound only
id=bc==2;
xu(id)=bounds(id,1)+exp(xt(id));
% upper bound only
id=bc==1;
xu(id)=bounds(id,2)-exp(xt(id));
end
