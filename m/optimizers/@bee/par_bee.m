classdef par_bee %< handle
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
        vargs={};
    end
    properties(SetAccess=protected)
        Objective
        best
        best_fval
        best_viol_strength
        best_fitness
        the_bees
        number_of_parameters
        finish_time
        optimizer='abc';
    end
    methods
        function obj=par_bee(Objective,x0,f0,lb,ub,options,varargin)
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
            %        RES = par_bee(@myfunc,...)
            %     In this case, F = myfunc(X) returns the scalar function value F of
            %     the MYFUNC function evaluated at X.
            %
            %     FUN can also be an anonymous function:
            %        RES = par_bee(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[],[0 0])
            %     returns X = [0;0].
            %
            %     FUN=inline('sum(x.^2)'); n=100;
            %     lb=-20*ones(n,1); ub=-lb; x0=lb+(ub-lb).*rand(n,1);
            %     optimpot=struct('MaxNodes',20,'MaxIter',1000,'MaxTime',60,...
            %     'MaxFunEvals',inf);
            %     RES=par_bee(@(x) FUN(x),x0,[],lb,ub,optimpot)
            
            % Reference: Inspired from Karaboga
            
            %   Copyright 2011 Junior Maih (junior.maih@gmail.com).
            %   $Revision: 8 $  $Date: 2011/08/13 11:23 $
            
            if nargin==0
                return
            end
            if nargin<5
                options=[];
            end
            if ~isfield(options,'nonlcon')||isempty(options.nonlcon)
                options.nonlcon=@(z)[lb-z;z-ub];
            else
                options.nonlcon=@(z)[options.nonlcon(z,varargin{:});lb-z;z-ub];
            end
            % the restrictions are satisfied when all elements in this
            % vector are negative or equal to 0
            % set these default options and change them if later on options is non-empty
            if ~isempty(options)
                all_prop=fieldnames(par_bee);
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
                obj.fitness_0=utils.optim.compute_fitness(obj.f0);
                obj=memorize_best_source(obj);
            end
            % Now find the peak
            obj=optimize_bees(obj);
        end
    end
end

function obj=optimize_bees(obj)
obj.food_number=round(.5*obj.MaxNodes);
obj.the_bees=struct('xx',{},'ff',{},'fitness',{},'trial',{},'violation_strength',{});
n0=size(obj.x0,2);
for i0=1:n0
    obj.the_bees(i0)=...
        struct('xx',obj.x0(:,i0),...
        'ff',obj.f0(i0),...
        'fitness',obj.fitness_0(i0),...
        'trial',0,...
        'violation_strength',obj.violation_strength_0(i0));
end
missing=obj.MaxNodes-n0;
% set and record the seed before we start drawing anything
s = RandStream.create('mt19937ar','seed',obj.rand_seed);
try
    RandStream.setGlobalStream(s);
catch %#ok<CTCH>
    RandStream.setDefaultStream(s); %#ok<SETRS>
end

[xx_,ff_,viol_strength_,funevals,fitness_,trial_]=...
    new_bees(obj.Objective,obj.lb,obj.ub,missing,obj.nonlcon,...
    obj.penalty,obj.vargs{:});
for i0=1:missing
    obj.the_bees(n0+i0)=...
        struct('xx',xx_(:,i0),...
        'ff',ff_(i0),...
        'fitness',fitness_(i0),...
        'trial',trial_(i0),...
        'violation_strength',viol_strength_(i0));
end
obj.funcCount=obj.funcCount+funevals;
obj=memorize_best_source(obj);

if ~obj.stopping_created
    utils.optim.manual_stopping;
else
    if ~exist('ManualStoppingFile.txt','file')
        utils.optim.manual_stopping;
    end
end
stopflag=utils.optim.check_convergence(obj);
while isempty(stopflag)
    obj.iterations=obj.iterations+1;
    % send employed bees
    %-------------------
    obj=send_bees(obj,'e');
    % send onlooker bees
    %-------------------
    obj=send_bees(obj,'o');
    % store best food source location
    %--------------------------------
    obj=memorize_best_source(obj);
    % send scout bees
    %----------------
    obj=send_bees(obj,'s');
    if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
        restart=1;
        fmin_iter=obj.best_fval;
        disperse=utils.optim.dispersion([obj.the_bees.xx],obj.lb,obj.ub);
        utils.optim.display_progress(restart,obj.iterations,obj.best_fval,fmin_iter,...
            disperse,obj.funcCount,obj.optimizer);
    end
    stopflag=utils.optim.check_convergence(obj);
end
obj.finish_time=clock;
end

function obj=send_bees(obj,type)
mut_func=@generate_mutant;
ZZ=[obj.the_bees.xx];
food_number=obj.food_number;
nparams=obj.number_of_parameters;
objfun=obj.Objective;
vargs=obj.vargs;
lb=obj.lb;
ub=obj.ub;
nonlcon=obj.nonlcon;
the_bees=obj.the_bees;
switch type
    case 's' % scout
        trials=[the_bees.trial];
        renew=find(trials>=obj.max_trials);
        [xx_,ff_,viol_strength_,funevals,fitness_,...
            trial_]=new_bees(objfun,lb,ub,numel(renew),nonlcon,obj.penalty,...
            vargs{:});
        iter=0;
        for ib=renew
            iter=iter+1;
            the_bees(ib)=struct('xx',xx_(:,iter),...
                'ff',ff_(iter),...
                'fitness',fitness_(iter),...
                'trial',trial_(iter),...
                'violation_strength',viol_strength_(iter));
        end
        obj.funcCount=obj.funcCount+funevals;
    case 'e' % employed bees
        parfor ii=1:food_number
            the_bees(ii)=mut_func(the_bees(ii),ii);
        end
        obj.funcCount=obj.funcCount+obj.food_number;
    case 'o' % onlooker bees
        fitness=[obj.the_bees.fitness];
        prob=(0.9.*fitness/max(fitness))+0.1;
        t=0;
        ii=1;
        while t<food_number
            if rand<prob(ii)
                t=t+1;
                the_bees(ii)=mut_func(the_bees(ii),ii);
                obj.funcCount=obj.funcCount+1;
            end
            ii=ii+1;
            if ii==food_number+1
                ii=1;
            end
        end
    otherwise
        error(['unknown type of bees ',type])
end
obj.the_bees=the_bees;
    function bb=generate_mutant(bb,ii)
        % generate a solution for index ii
        % a randomly chosen solution different from ii is used
        % for producing a mutant
        while 1
            donor_id=min(food_number,fix(rand*food_number)+1);
            if donor_id~=ii
                break
            end
        end
        mutant=bb.xx;
        % %     pick the parameters to change in the new solution
        change=min(fix(rand*nparams)+1,nparams);
        mutant(change)=mutant(change)+(mutant(change)-...
            ZZ(change,donor_id))*2.*(rand(numel(change),1)-.5);
        mutant(change)=utils.optim.recenter(mutant(change),lb(change),ub(change));
        f_mut=objfun(mutant,vargs{:});
        viol=nonlcon(mutant);
        viol=viol(viol>0);
        viol_strength=sum(viol);
        fit_mut=utils.optim.compute_fitness(f_mut);
        bb=deb_selection(bb,viol_strength,mutant,fit_mut,f_mut);
    end
end

function obj=memorize_best_source(obj)
if isempty(obj.the_bees)
    obj.best_fitness=obj.fitness_0(1);
    obj.best_fval=obj.f0(1);
    obj.best=obj.x0(:,1);
    obj.best_viol_strength=obj.violation_strength_0(1);
else
    obj=sort_bees(obj);
    tmp=deb_selection(obj.the_bees(1),obj.best_viol_strength,obj.best,obj.best_fitness,obj.best_fval);
    if isempty(obj.best_fval)||tmp.ff<obj.best_fval
        obj.best_fval=tmp.ff;
        obj.best=tmp.xx;
        obj.best_viol_strength=tmp.violation_strength;
    end
end
end

function obj=sort_bees(obj)
VS=unique([obj.the_bees.violation_strength]);
bb=struct('xx',{},'ff',{},'fitness',{},'trial',{},'violation_strength',{});
for iloc=1:numel(VS)
    this=[obj.the_bees.violation_strength]==VS(iloc);
    nthis=sum(this); % this is logical and so the sum has to be taken
    [~,tag]=sort([obj.the_bees(this).ff]);
    bb(end+(1:nthis))=obj.the_bees(tag);
end
obj.the_bees=bb;
end

function [x,f,viol_strength,funevals,fit,trial]=new_bees(objective,lb,ub,n,nonlcon,penalty,...
    varargin)
[x,f,viol_strength,funevals]=utils.optim.generate_candidates(objective,lb,ub,n,nonlcon,penalty,varargin{:});
trial=zeros(1,n);
fit=utils.optim.compute_fitness(f);
for ii=1:n
    v=viol_strength{ii};
    v=v(v>0);
    viol_strength{ii}=sum(v);
end
viol_strength=cell2mat(viol_strength);
end


function bb=deb_selection(bb,viol_strength,mutant,fit_mut,f_mut)
% we apply a greedy selection between ii and the mutant after controling
% for feasibility
if (viol_strength<bb.violation_strength)||...
        (viol_strength==bb.violation_strength && fit_mut>bb.fitness)
    bb.fitness=fit_mut;
    bb.xx=mutant;
    bb.ff=f_mut;
    bb.trial=0;
    bb.violation_strength=viol_strength;
else
    bb.trial=bb.trial+1;
end
end

