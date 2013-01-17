classdef bee %< handle
    properties
        stopping_created=false;
        start_time
        MaxNodes=20;
        lb
        ub
        x0
        f0
        violation_strength_0
        fitness_0
        iterations=0;
        funcCount=0;
        MaxIter=1000;
        MaxTime=3600;
        MaxFunEvals=inf;
        rand_seed=100*sum(clock);
        penalty=1e+8;
        verbose=10;
        % optimizer-specific properties
        max_trials=100;
        nonlcon
    end
    %     properties(Constant, Hidden = true)
    %     end
    properties (SetAccess = private, Hidden = true)
        food_number
        trial
        fitness
        vargs={};
    end
    properties(SetAccess=protected)
        Objective
        best
        best_fval
        best_viol_strength
        best_fitness
        xx
        ff
        number_of_parameters
        finish_time
        optimizer='abc';
        violation_strength
    end
    methods
        function obj=bee(Objective,x0,f0,lb,ub,options,varargin)
            %  BEE attempts to find the global minimum of a constrained function of
            %  several variables.
            %   BEE attempts to solve problems of the form:
            %    min F(X)  subject to:  LB <= X <= UB   (bounds)
            %     X
            %
            %   RES = BEE constructs an object with the default optimization parameters.
            %   RES has fields that are intuitive to understand. 'lb' (lower bound),
            %   'ub' (upper bound),'x0',(vector of initial values),'f0'(function value
            %   at x0), 'vargs'(additional arguments of FUN),'penalty' (threshold
            %   function value beyond which the parameter draws are
            %   discarded),'Objective' (name of the objective function), 'best'(best
            %   parameter vector), 'best_fval'(best function value), 'xx'(parameter
            %   vectors in the colony),'ff'(function values at xx)
            %   'MaxIter','MaxTime','MaxNodes
            %
            %   RES = BEE(FUN,X0,[],LB,UB) starts at X0 and finds a minimum X to the
            %   function FUN, subject to the bound constraints LB and UB. FUN accepts
            %   input X and returns a scalar function value F evaluated at X. X0 may be
            %   a scalar or a vector.
            %
            %   RES = BEE(FUN,X0,[],LB,UB,OPTIONS) optimizes the function FUN under the
            %   optimization options set under the structure OPTIONS. The fields of
            %   this structure could be all or any of the following:
            %       - 'MaxNodes': this the number of different elements in the group
            %       that will share information with each other in order to find better
            %       solutions. The default is 20
            %       - 'MaxIter': the maximum number of iterations. The default is 1000
            %       - 'MaxTime': The time budget in seconds. The default is 3600
            %       - 'MaxFunEvals': the maximum number of function evaluations. The
            %       default is inf
            %       - 'rand_seed': the seed number for the random draws
            %
            %   Optimization stops when one of the following happens:
            %   1- the number of iterations exceeds MaxIter
            %   2- the number of function counts exceeds MaxFunEvals
            %   3- the time elapsed exceeds MaxTime
            %   4- the user write anything in and saves the automatically generated
            %   file called "ManualStopping.txt"
            %
            %   Examples
            %     FUN can be specified using @:
            %        RES = bee(@myfunc,...)
            %     In this case, F = myfunc(X) returns the scalar function value F of
            %     the MYFUNC function evaluated at X.
            %
            %     FUN can also be an anonymous function:
            %        RES = bee(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[],[0 0])
            %     returns X = [0;0].
            %
            %     FUN=inline('sum(x.^2)'); n=100;
            %     lb=-20*ones(n,1); ub=-lb; x0=lb+(ub-lb).*rand(n,1);
            %     optimpot=struct('MaxNodes',20,'MaxIter',1000,'MaxTime',60,...
            %     'MaxFunEvals',inf);
            %     RES=bee(@(x) FUN(x),x0,[],lb,ub,optimpot)
            
            % Reference: Inspired from Karaboga
            
            %   Copyright 2011 Junior Maih (junior.maih@gmail.com).
            %   $Revision: 8 $  $Date: 2011/08/13 11:23 $
            
            if nargin==0
                return
            end
            if nargin<5
                options=[];
            end
            if isempty(options.nonlcon)
                options.nonlcon=@(z)0;
            end
            options.nonlcon=@(z)[options.nonlcon(z);lb-z;z-ub];
            % the restrictions are satisfied when all elements in this
            % vector are negative or equal to 0
            % set these default options and change them if later on options is non-empty
            if ~isempty(options)
                all_prop=fieldnames(bee);
                options_prop=fieldnames(options);
                nargs=numel(options_prop);
                for ii=1:nargs
                    prop=options_prop{ii};
                    if isempty(find(strcmp(prop,all_prop),1))
                        error([mfilename,':: ',prop,...
                            ' is not a valid property of ',...
                            mfilename,' objects'])
                    end
                    obj.(prop)=options.(prop);
                end
            end
            obj.lb=lb;
            obj.ub=ub;
            obj.x0=x0;
            obj.f0=f0;
            if isempty(obj.start_time)
                obj.start_time=clock;
            else
                if numel(obj.start_time)~=6
                    error([mfilename,':: wrong entry for start_time (should be same format as clock)'])
                end
            end
            obj.vargs=varargin;
            obj.number_of_parameters=size(obj.lb,1);
            obj.Objective=fcnchk(Objective,length(obj.vargs));
            n0=size(obj.x0,2);
            if n0
                n0=min(n0,obj.MaxNodes);
                obj.x0=obj.x0(:,1:n0);
                obj.violation_strength_0=zeros(1,n0);
                obj.f0=nan(1,1:n0);
                for ii=1:n0
                    obj.f0(ii)=obj.Objective(obj.x0(:,ii),obj.vargs{:});
                    tmp=obj.nonlcon(obj.x0(:,ii));
                    tmp=tmp(tmp>0);
                    obj.violation_strength_0(ii)=sum(tmp);
                end
                obj.funcCount=obj.funcCount+n0;
                obj.fitness_0=compute_fitness(obj.f0);
                obj=memorize_best_source(obj);
            end
            % Now find the peak
            obj=optimize_bees(obj);
        end
    end
end

function obj=optimize_bees(obj)
obj.food_number=round(.5*obj.MaxNodes);
obj.xx=nan(obj.number_of_parameters,obj.MaxNodes);
obj.ff=nan(1,obj.MaxNodes);
obj.fitness=nan(1,obj.MaxNodes);
obj.trial=nan(1,obj.MaxNodes);
obj.violation_strength=nan(1,obj.MaxNodes);
n0=size(obj.x0,2);
if n0
    obj.fitness(1:n0)=obj.fitness_0(1:n0);
    obj.trial(1:n0)=0;
    obj.ff(1:n0)=obj.f0;
    obj.xx(:,1:n0)=obj.x0(:,1:n0);
    obj.violation_strength(1:n0)=obj.violation_strength_0;
end
missing=obj.MaxNodes-n0;
% set and record the seed before we start drawing anything
s = RandStream.create('mt19937ar','seed',obj.rand_seed);
try
    RandStream.setGlobalStream(s);
catch %#ok<CTCH>
    RandStream.setDefaultStream(s); %#ok<SETRS>
end

[obj.xx(:,n0+1:end),obj.ff(n0+1:end),...
    obj.violation_strength(n0+1:end),funevals,...
    obj.fitness(n0+1:end),obj.trial(n0+1:end)]=...
    new_bees(obj.Objective,obj.lb,obj.ub,missing,obj.nonlcon,...
    obj.penalty,obj.vargs{:});
obj.funcCount=obj.funcCount+funevals;
obj=memorize_best_source(obj);

if ~obj.stopping_created
    manual_stopping;
else
    if ~exist('ManualStoppingFile.txt','file')
        manual_stopping;
    end
end
stopflag=check_convergence(obj);
while isempty(stopflag)
    obj.iterations=obj.iterations+1;
    obj=send_employed_bees(obj);
    obj=send_onlooker_bees(obj);
    obj=memorize_best_source(obj);
    obj=send_scout_bees(obj);
    if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
        restart=1;
        fmin_iter=obj.best_fval;
        disperse=dispersion(obj.xx,obj.lb,obj.ub);
        display_progress(restart,obj.iterations,obj.best_fval,fmin_iter,...
            disperse,obj.funcCount,obj.optimizer);
    end
    stopflag=check_convergence(obj);
end
obj.finish_time=clock;
end

function obj=send_scout_bees(obj)
renew=find(obj.trial>=obj.max_trials);
[obj.xx(:,renew),obj.ff(renew),obj.violation_strength(renew),funevals,obj.fitness(renew),...
    obj.trial(renew)]=new_bees(obj.Objective,obj.lb,...
    obj.ub,numel(renew),obj.nonlcon,obj.penalty,obj.vargs{:});
obj.funcCount=obj.funcCount+funevals;
end

function obj=send_employed_bees(obj)
for ii=1:obj.food_number
    obj=generate_mutant(obj,ii);
end
end

function obj=send_onlooker_bees(obj)
prob=(0.9.*obj.fitness./max(obj.fitness))+0.1;
t=0;
ii=1;
while t<obj.food_number
    if rand<prob(ii)
        t=t+1;
        obj=generate_mutant(obj,ii);
    end
    ii=ii+1;
    if ii==obj.food_number+1
        ii=1;
    end
end
end

function obj=memorize_best_source(obj)
if isempty(obj.ff)
    obj.best_fitness=obj.fitness_0(1);
    obj.best_fval=obj.f0(1);
    obj.best=obj.x0(:,1);
    obj.best_viol_strength=obj.violation_strength_0(1);
else
    obj=sort_bees(obj);
    tmp=deb_selection(obj,1,obj.best_viol_strength,obj.best,obj.best_fitness,obj.best_fval);
    if isempty(obj.best_fval)||obj.ff(1)<obj.best_fval
        obj.best_fval=tmp.ff(1);
        obj.best=tmp.xx(:,1);
        obj.best_viol_strength=obj.violation_strength(1);
    end
end
end

function obj=sort_bees(obj)
VS=unique(obj.violation_strength);
FF=nan(1,obj.MaxNodes);
XX=nan(obj.number_of_parameters,obj.MaxNodes);
FIT=nan(1,obj.MaxNodes);
TRIALS=nan(1,obj.MaxNodes);
VIOLS=nan(1,obj.MaxNodes);
iter0=0;
for iloc=1:numel(VS)
    this=obj.violation_strength==VS(iloc);
    nthis=numel(this);
    [FF(iter0+(1:nthis)),XX(:,iter0+(1:nthis)),FIT(iter0+(1:nthis)),TRIALS(iter0+(1:nthis)),VIOLS(iter0+(1:nthis))]=...
        resort(obj.ff(this),obj.xx(:,this),obj.fitness(this),obj.trial(this),obj.violation_strength(this));
    iter0=iter0+nthis;
end
[obj.ff,obj.xx,obj.fitness,obj.trial,obj.violation_strength]=deal(FF,XX,FIT,TRIALS,VIOLS);
end

function [x,f,viol_strength,funevals,fit,trial]=new_bees(objective,lb,ub,n,nonlcon,penalty,...
    varargin)
[x,f,viol_strength,funevals]=generate_candidates(objective,lb,ub,n,nonlcon,penalty,varargin{:});
trial=zeros(1,n);
fit=compute_fitness(f);
for ii=1:n
    v=viol_strength{ii};
    v=v(v>0);
    viol_strength{ii}=sum(v);
end
viol_strength=cell2mat(viol_strength);
end

function obj=deb_selection(obj,ii,viol_strength,mutant,fit_mut,f_mut)
% we apply a greedy selection between ii and the mutant after controling
% for feasibility
if (viol_strength<obj.violation_strength(ii))||...
        (viol_strength==obj.violation_strength(ii) && fit_mut>obj.fitness(ii))
    obj.fitness(ii)=fit_mut;
    obj.xx(:,ii)=mutant;
    obj.ff(:,ii)=f_mut;
    obj.trial(ii)=0;
    obj.violation_strength(ii)=viol_strength;
else
    obj.trial(ii)=obj.trial(ii)+1;
end
end

function obj=generate_mutant(obj,ii)
% generate a solution for index ii
% a randomly chosen solution different from ii is used
% for producing a mutant
while 1
    donor_id=min(obj.food_number,fix(rand*obj.food_number)+1);
    if donor_id~=ii
        break
    end
end
mutant=obj.xx(:,ii);

% %     pick the parameters to change in the new solution
change=min(fix(rand*obj.number_of_parameters)+1,obj.number_of_parameters);
mutant(change)=mutant(change)+(mutant(change)-...
    obj.xx(change,donor_id))*2.*(rand(numel(change),1)-.5);
mutant(change)=recenter(mutant(change),obj.lb(change),obj.ub(change));
f_mut=obj.Objective(mutant,obj.vargs{:});
viol=obj.nonlcon(mutant);
viol=viol(viol>0);
viol_strength=sum(viol);
obj.funcCount=obj.funcCount+1;
fit_mut=compute_fitness(f_mut);
obj=deb_selection(obj,ii,viol_strength,mutant,fit_mut,f_mut);
end
