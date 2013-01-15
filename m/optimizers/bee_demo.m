classdef bee_demo < handle
    properties
        stopping_created=false;
        start_time
        colony_size=20;
        lb
        ub
        x0
        f0
        max_iter=1000;
        max_time=3600;
        max_fcount=inf;
        rand_seed=100*sum(clock);
        penalty=1e+8;
        display=10;
        % optimizer-specific properties
        max_trials=100;
        fbest
        xbest
        speed=.25;
        griddata
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
        xx
        ff
        number_of_parameters
        finish_time
        iter=0;
        fcount=0;
        optimizer='abc';
    end
    methods(Access=private)
        function obj=optimize(obj)
            obj.food_number=round(.5*obj.colony_size);
            obj.xx=nan(obj.number_of_parameters,obj.colony_size);
            obj.ff=nan(1,obj.colony_size);
            obj.fitness=nan(1,obj.colony_size);
            obj.trial=nan(1,obj.colony_size);
            n0=size(obj.x0,2);
            obj.fitness(1:n0)=compute_fitness(obj.f0(1:n0));
            obj.trial(1:n0)=0;
            obj.ff(1:n0)=obj.f0;
            obj.xx(:,1:n0)=obj.x0(:,1:n0);
            missing=obj.colony_size-n0;
            % set and record the seed before we start drawing anything
            s = RandStream.create('mt19937ar','seed',obj.rand_seed);
            RandStream.setDefaultStream(s);
            
            [obj.xx(:,n0+1:end),obj.ff(n0+1:end),funevals,...
                obj.fitness(n0+1:end),obj.trial(n0+1:end)]=...
                new_bees(obj.Objective,obj.lb,obj.ub,missing,...
                obj.penalty,obj.vargs{:});
            obj.fcount=obj.fcount+funevals;
            obj.memorize_best_source;
            
            if ~obj.stopping_created
                manual_stopping;
            else
                if ~exist('ManualStoppingFile.txt','file')
                    manual_stopping;
                end
            end
            hh=[];
            stopflag=check_convergence(obj);
            while isempty(stopflag)
                obj.iter=obj.iter+1;
                obj.send_employed_bees;
                obj.send_onlooker_bees;
                obj.memorize_best_source;
                obj.send_scout_bees;
                if rem(obj.iter,obj.display)==0 || obj.iter==1
                    restart=1;
                    fmin_iter=obj.best_fval;
                    disperse=dispersion(obj.xx,obj.lb,obj.ub);
                    display_progress(restart,obj.iter,obj.best_fval,fmin_iter,...
                        disperse,obj.fcount,obj.optimizer);
                end
                stopflag=check_convergence(obj);
                hh=vizualize_progress(hh,mfilename,obj.xx,obj.ff,obj.xbest,obj.fbest,obj.griddata,obj.lb,obj.ub);
                pause(obj.speed)
            end
            obj.finish_time=clock;
        end
        
        function obj=send_scout_bees(obj)
            renew=find(obj.trial>=obj.max_trials);
            test=true;
            if test
                [obj.xx(:,renew),obj.ff(renew),funevals,obj.fitness(renew),...
                    obj.trial(renew)]=new_bees4_renewal(obj.Objective,...
                    obj.best,obj.lb,obj.ub,numel(renew),obj.vargs{:});
            else
                [obj.xx(:,renew),obj.ff(renew),funevals,obj.fitness(renew),...
                    obj.trial(renew)]=new_bees(obj.Objective,obj.lb,...
                    obj.ub,numel(renew),obj.penalty,obj.vargs{:});
            end
            obj.fcount=obj.fcount+funevals;
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
            [obj.ff,obj.xx,obj.fitness,obj.trial]=...
                resort(obj.ff,obj.xx,obj.fitness,obj.trial);
%             [obj.ff,obj.xx,obj.fitness]=...
%                 resort(obj.ff,obj.xx,obj.fitness);
            if isempty(obj.best_fval)||obj.ff(1)<obj.best_fval
                obj.best_fval=obj.ff(1);
                obj.best=obj.xx(:,1);
            end
        end
    end
    methods
        function obj=bee_demo(Objective,x0,f0,lb,ub,options,varargin)
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
            %   'max_iter','max_time','colony_size
            %
            %   RES = BEE(FUN,X0,[],LB,UB) starts at X0 and finds a minimum X to the
            %   function FUN, subject to the bound constraints LB and UB. FUN accepts
            %   input X and returns a scalar function value F evaluated at X. X0 may be
            %   a scalar or a vector.
            %
            %   RES = BEE(FUN,X0,[],LB,UB,OPTIONS) optimizes the function FUN under the
            %   optimization options set under the structure OPTIONS. The fields of
            %   this structure could be all or any of the following:
            %       - 'colony_size': this the number of different elements in the group
            %       that will share information with each other in order to find better
            %       solutions. The default is 20
            %       - 'max_iter': the maximum number of iterations. The default is 1000
            %       - 'max_time': The time budget in seconds. The default is 3600
            %       - 'max_fcount': the maximum number of function evaluations. The
            %       default is inf
            %       - 'rand_seed': the seed number for the random draws
            %
            %   Optimization stops when one of the following happens:
            %   1- the number of iterations exceeds max_iter
            %   2- the number of function counts exceeds max_fcount
            %   3- the time elapsed exceeds max_time
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
            %     optimpot=struct('colony_size',20,'max_iter',1000,'max_time',60,...
            %     'max_fcount',inf);
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
            % set these default options and change them if later on options is non-empty
            if ~isempty(options)
                all_prop=fieldnames(bee_demo);
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
                n0=min(n0,obj.colony_size);
                obj.x0=obj.x0(:,1:n0);
                if isempty(obj.f0)
                    obj.f0=nan(1,1:n0);
                    for ii=1:n0
                        obj.f0(ii)=obj.Objective(obj.x0(:,ii),obj.vargs{:});
                    end
                    obj.fcount=obj.fcount+n0;
                else
                    obj.f0=obj.f0(1:n0);
                end
                [obj.f0,obj.x0]=resort(obj.f0,obj.x0);
                obj.best_fval=obj.f0(1);
                obj.best=obj.x0(:,1);
            end
            
            % Now find the peak
            obj.optimize;
        end
    end
end

function [x,f,funevals,fit,trial]=new_bees(objective,lb,ub,n,penalty,...
    varargin)
[x,f,funevals]=generate_candidates(objective,lb,ub,n,penalty,varargin{:});
trial=zeros(1,n);
fit=compute_fitness(f);
end

function [x,f,funevals,fit,trial]=new_bees4_renewal(objective,xbest,lb,ub,...
    n,varargin)
npar=size(lb,1);
x=bsxfun(@plus,lb,bsxfun(@times,ub-lb,rand(npar,n)));
x=x+rand(1,1).*(xbest(:,ones(1,n))-x);
x=recenter(x,lb,ub);
f=nan(1,n);
for ii=1:n
    f(ii)=objective(x(:,ii),varargin{:});
end
trial=zeros(1,n);
fit=compute_fitness(f);
funevals=n;
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
% pick the parameters to change in the new solution
change=min(fix(rand*obj.number_of_parameters)+1,obj.number_of_parameters);
mutant(change)=mutant(change)+(mutant(change)-...
    obj.xx(change,donor_id))*2*(rand-.5);
mutant=recenter(mutant,obj.lb,obj.ub);
f_mut=obj.Objective(mutant,obj.vargs{:});
obj.fcount=obj.fcount+1;
fit_mut=compute_fitness(f_mut);
% we apply a greedy selection between ii and the
% mutant
if fit_mut>obj.fitness(ii)
    obj.fitness(ii)=fit_mut;
    obj.xx(:,ii)=mutant;
    obj.ff(:,ii)=f_mut;
    obj.trial(ii)=0;
else
    obj.trial(ii)=obj.trial(ii)+1;
end
end
