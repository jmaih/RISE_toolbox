function [xbest,fbest,H,issue]=estimation_engine(PROBLEM,estim_blocks)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

issue='';

mydefaults=the_defaults();
    
if nargin==0
    
    if nargout
        
        xbest=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if nargin<2
    
    estim_blocks=[];
    
end

opt=disp_defaults(mydefaults);

OtherProblemFields={'Aineq','bineq','Aeq','beq','nonlcon'};

if isempty(PROBLEM.solver)
    
    PROBLEM.solver=opt.optimizer;
    
end

for ii=1:numel(OtherProblemFields)
    
    if ~isfield(PROBLEM,OtherProblemFields{ii})
        
        PROBLEM.(OtherProblemFields{ii})=[];
        
    end
    
end

if isnumeric(PROBLEM.solver)
    
    error('optimizer should be a string or a function handle')
    
end

optimizer=PROBLEM.solver;

fh=PROBLEM.objective;

x0=PROBLEM.x0;

lb=PROBLEM.lb;

ub=PROBLEM.ub;

optim_opt=PROBLEM.options;

options=opt.optimset;

options=utils.miscellaneous.setfield(options,optim_opt);

block_flag=~isempty(estim_blocks);

if block_flag
    
    options.blocks=estim_blocks;
    
end

vargs={};

if iscell(optimizer)
    
    vargs=optimizer(2:end);
    
    optimizer=optimizer{1};
    
end

% rename the optimizers if necessary
if isa(optimizer,'function_handle')
    
    tmp=func2str(optimizer);
    
else
    
    tmp=optimizer;
    
end

if ~strcmp(tmp,'fmincon')
    
    options.nonlcon=PROBLEM.nonlcon;
    
end

if ismember(tmp,{'fminunc','fminsearch'})
    
    optimizer=str2func(['@',tmp,'_bnd']);
    
end

if ischar(optimizer)
    
    optimizer=str2func(optimizer);
    
end

% hessian
%---------
d=numel(x0);

H=nan(d);

if block_flag
    
    options.optimizer=optimizer;
    
    [xfinal,ffinal,exitflag]=blockwise_optimization(fh,x0,lb,ub,options,vargs{:});

else
    
    if strcmp(tmp,'fmincon')
        
        [xfinal,ffinal,exitflag,~,~,~,H(:,:,1)]=optimizer(fh,x0,[],[],[],[],lb,ub,PROBLEM.nonlcon,options,vargs{:});
    
    elseif strcmp(tmp,'fminunc')
        
        [xfinal,ffinal,exitflag,~,~,H(:,:,1)]=optimizer(fh,x0,lb,ub,options,vargs{:});
    
    elseif strcmp(tmp,'fminsearch')
        
        [xfinal,ffinal,exitflag]=optimizer(fh,x0,lb,ub,options,vargs{:});
        
    else
        
        nout=nargout(optimizer);
        
        if nout==3
            
            [xfinal,ffinal,exitflag]=optimizer(fh,x0,lb,ub,options,vargs{:});
            
        elseif nout>=4
            
            [xfinal,ffinal,exitflag,H]=optimizer(fh,x0,lb,ub,options,vargs{:});
        
        else
            
            error('optimizer should have at least 3 outputs')
            
        end
        
    end
    
end

% select the overall best
xbest=xfinal;

fbest=ffinal;

TranslateOptimization(exitflag);

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

function [x,fval,exitflag,output]=fminsearch_bnd(fun,x0,lb,ub,options,varargin) %#ok<DEFNU>

nonlcon=options.nonlcon;

options=rmfield(options,'nonlcon');

bc=bound_class([lb,ub]);

xt=transform_parameters(x0,bc,[lb,ub]);

[xt,fval,exitflag,output]=fminsearch(@wrapper,xt,options,fun,bc,lb,ub,varargin{:});

x=untransform_parameters(xt,bc,[lb,ub]);

end

function [x,fval,exitflag,output,grad,hess]=fminunc_bnd(fun,x0,lb,ub,options,varargin) %#ok<DEFNU>

nonlcon=options.nonlcon;

options=rmfield(options,'nonlcon');

bc=bound_class([lb,ub]);

xt=transform_parameters(x0,bc,[lb,ub]);

[xt,fval,exitflag,output,grad,hess]=fminunc(@wrapper,xt,options,fun,bc,lb,ub,varargin{:});

x=untransform_parameters(xt,bc,[lb,ub]);

end

function ff=wrapper(xt,fun,bc,lb,ub,varargin)

xu=untransform_parameters(xt,bc,[lb,ub]);

ff=fun(xu,varargin{:});

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
% 'Algorithm','interior-point',...% [ active-set | interior-point | levenberg-marquardt | sqp | ...
% I observed the following problem with the interior point algorithm (IP):
% I estimated a model using my abc routines and then wanted to sharpen the
% results using fmincon. fmincon with the IP could not even replicate the
% function at the starting point (my guess is that this is because some
% parameters were lying up against their boundaries). More exactly, the
% function was given the value of the penalty. When I turned off the
% interior point, everything was fine. Maybe all this is consistent with
% the definition of interior point. Further testing with wrapper function
% fminsearch_bnd also gave the wrong starting value function. This points
% in the direction of a penalty arising with the parameters not being in
% the interior.
% I probably should revisit the way I impose constraints in abc: this
% should be obsolete by now

function d=the_defaults()

opt.optimset=optimset('Display','iter',...%[ off | iter | iter-detailed | notify | notify-detailed | final | final-detailed ]
    'MaxFunEvals',inf,...% [ positive scalar ]
    'MaxIter',1000,...%: [ positive scalar ]
    'TolFun',sqrt(eps),...%: [ positive scalar ]
    'TolX',sqrt(eps),...%: [ positive scalar ]
    'MaxNodes',20,...%: [ positive scalar | {1000*numberOfVariables} ]
    'UseParallel','never',...%: [ always | {never} ]
    'MaxTime',3600);%: [ positive scalar | {7200} ]

opt.optimizer='fmincon'; % default is fmincon

d={
    'optimset(sr)',opt.optimset,@(x)isstruct(x),...
    'optimset must be a structure'
    
    'optimizer',opt.optimizer,...
    @(x)ischar(x)||iscell(x)||isa(x,'function_handle'),...
    'optimizer must be a char, a cell or a function handle'
    
    };

end

