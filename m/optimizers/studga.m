classdef studga < handle
    
    properties
        stopping_created=false;
        start_time
        colony_size=50;
        lb
        ub
        x0
        f0
        max_iter=1000;
        max_time=3600;
        max_fcount=inf;
        rand_seed=100*sum(clock);
        penalty=1e+8;
        verbose=10;
        % optimizer-specific options
        crossover_type = 1; % crossover type: 1 = single point, 2 = two point, 3 = uniform
        crossover_probability = 1; % crossover probability
        mutation_probability = 0.01; % initial mutation probability
        elitism=2;
		restrictions
    end
    %     properties(Constant, Hidden = true)
    %     end
    properties (SetAccess = private, Hidden = true)
        vargs={};
    end
    properties(SetAccess=protected)
        Objective
        best
        best_fit
        xx
        ff
        number_of_parameters
        finish_time
        iter=0;
        fcount=0;
        optimizer='studga';
    end
    methods(Access=private)
        function obj=optimize(obj)
            obj.xx=nan(obj.number_of_parameters,obj.colony_size);
            obj.ff=nan(1,obj.colony_size);
            n0=size(obj.x0,2);
            if n0
                obj.ff(1:n0)=obj.f0;
                obj.xx(:,1:n0)=obj.x0(:,1:n0);
            end
            missing=obj.colony_size-n0;
            % set and record the seed before we start drawing anything
            s = RandStream.create('mt19937ar','seed',obj.rand_seed);
            RandStream.setDefaultStream(s);
            
            [obj.xx(:,n0+1:end),obj.ff(n0+1:end),funevals]=...
                generate_candidates(obj.Objective,obj.lb,obj.ub,missing,...
                obj.restrictions,obj.penalty,obj.vargs{:});
            [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
            % % %             obj.probabilities=1/obj.colony_size*ones(1,obj.colony_size);
            obj.fcount=obj.fcount+funevals;
            obj.memorize_best_solution;
            
            %==================
            if isnumeric(obj.crossover_type)
                obj.crossover_type=round(obj.crossover_type);
                if obj.crossover_type<=0 || obj.crossover_type>obj.number_of_parameters
                    error([mfilename,':: crossover_type be between 1 and ',int2str(obj.number_of_parameters)])
                end
            end
            %==================
            
            if ~obj.stopping_created
                manual_stopping;
            else
                if ~exist('ManualStoppingFile.txt','file')
                    manual_stopping;
                end
            end
            stopflag=check_convergence(obj);
            while isempty(stopflag)
                obj.iter=obj.iter+1;
                obj.crossover_selection;
                
                % Mutation
                obj.mutation;
                
                % Make sure the population does not have duplicates.
                obj.xx=clear_duplicates(obj.xx,obj.lb,obj.ub);
                
                % Make sure each individual is legal.
                obj.xx=recenter(obj.xx,obj.lb,obj.ub);
                
                % Calculate cost
                for k = 1 : obj.colony_size
                    obj.ff(k) = obj.Objective(obj.xx(:,k),obj.vargs{:});
                end
                obj.fcount=obj.fcount+obj.colony_size;
                
                % Sort from best to worst
                [obj.ff,obj.xx] = resort(obj.ff,obj.xx);
                
                obj.memorize_best_solution;
                if rem(obj.iter,obj.verbose)==0 || obj.iter==1
                    restart=1;
                    fmin_iter=obj.best_fit;
                    disperse=dispersion(obj.xx,obj.lb,obj.ub);
                    display_progress(restart,obj.iter,obj.best_fit,fmin_iter,...
                        disperse,obj.fcount,obj.optimizer);
                end
                stopflag=check_convergence(obj);
            end
            obj.finish_time=clock;
        end
        
        function obj=crossover_selection(obj)
            InverseCost = compute_fitness(obj.ff);
            %             InverseCost = 1./obj.ff;
            children=nan(obj.number_of_parameters,obj.colony_size-obj.elitism);
            for k = obj.elitism+1:2:obj.colony_size % begin selection/crossover loop
                % Select the most fit individual (the stud) as the first parent
                mate =[1,nan];
                % Select another parent to mate with the stud and create two children - roulette wheel selection
                % Make sure the stud is not selected as the second parent
                Random_Cost = 0;
                while InverseCost(1) >= Random_Cost
                    Random_Cost = rand*sum(InverseCost);
                end
                Select_Cost = InverseCost(1);
                Select_index = 1;
                while Select_Cost < Random_Cost
                    Select_index = Select_index + 1;
                    if Select_index >= obj.colony_size
                        break;
                    end
                    Select_Cost = Select_Cost + InverseCost(Select_index);
                end
                mate(2) = Select_index;
                Parents = obj.xx(:,mate);
                % Crossover
                % clone the parents and change children only if necessary
                children(:,k-obj.elitism) = Parents(:,1);
                children(:,k-obj.elitism+1) = Parents(:,2);
                if isnumeric(obj.crossover_type)||...
                        (ischar(obj.crossover_type) && strcmpi(obj.crossover_type,'variable'))
                    if isnumeric(obj.crossover_type)
                        npoints=obj.crossover_type;
                    else
                        npoints=ceil(rand*obj.number_of_parameters);
                    end
                    if obj.crossover_probability > rand
                        Xover_Pts = ceil(rand(1,npoints)*obj.number_of_parameters);
                        Intervals=[unique([1,Xover_Pts]),obj.number_of_parameters+1];
                        odd=true;
                        for pp=1:numel(Intervals)-1
                            ss=Intervals(pp):Intervals(pp+1)-1;
                            if odd
                                children(ss,k-obj.elitism) = Parents(ss,1);
                                children(ss,k-obj.elitism+1)= Parents(ss,2);
                            else
                                children(ss,k-obj.elitism) = Parents(ss,2);
                                children(ss,k-obj.elitism+1)= Parents(ss,1);
                            end
                            odd=~odd;
                        end
                    end
                elseif ischar(obj.crossover_type) && strcmpi(obj.crossover_type,'uniform')
                    for ii = 1 : obj.number_of_parameters
                        if obj.crossover_probability > rand
                            children(ii,k-obj.elitism) = Parents(ii,1);
                            children(ii,k-obj.elitism+1) = Parents(ii,2);
                        end
                    end
                else
                    error([mfilename,':: unknown type of crossover ',obj.crossover_type])
                end
            end % end selection/crossover loop
            % Replace the non-elite population members with the new children
            for k = obj.elitism+1:2:obj.colony_size
                obj.xx(:,k) = children(:,k-obj.elitism);
                obj.xx(:,k+1) = children(:,k-obj.elitism+1);
            end
        end
        
        function obj=mutation(obj)
            for individual = obj.elitism+1:obj.colony_size % Don't allow the elites to be mutated
                for parnum = 1 : obj.number_of_parameters
                    if obj.mutation_probability > rand
                        obj.xx(parnum,individual) = obj.lb(parnum) + (obj.ub(parnum)-obj.lb(parnum)) * rand;
                    end
                end
            end
        end
        
        function obj=memorize_best_solution(obj)
            [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
            if isempty(obj.best_fit)||obj.ff(1)<obj.best_fit
                obj.best_fit=obj.ff(1);
                obj.best=obj.xx(:,1);
            end
        end
    end
    methods
        function obj=studga(Objective,x0,f0,lb,ub,options,varargin)
            %  STUDGA attempts to find the global minimum of a constrained function of
            %  several variables
            %   STUDGA attempts to solve problems of the form:
            %    min F(X)  subject to:  LB <= X <= UB   (bounds)
            %     X
            %
            %   RES = STUDGA constructs an object with the default optimization parameters.
            %   RES has fields that are intuitive to understand. 'lb' (lower bound),
            %   'ub' (upper bound),'x0',(vector of initial values),'f0'(function value
            %   at x0), 'vargs'(additional arguments of FUN),'penalty' (threshold
            %   function value beyond which the parameter draws are
            %   discarded),'Objective' (name of the objective function), 'best'(best
            %   parameter vector), 'best_fit'(best function value), 'xx'(parameter
            %   vectors in the colony),'ff'(function values at xx)
            %   'max_iter','max_time','colony_size
            %
            %   RES = STUDGA(FUN,X0,LB,UB) starts at X0 and finds a minimum X to the
            %   function FUN, subject to the bound constraints LB and UB. FUN accepts
            %   input X and returns a scalar function value F evaluated at X. X0 may be
            %   a scalar or a vector.
            %
            %   RES = STUDGA(FUN,X0,LB,UB,OPTIONS) optimizes the function FUN under the
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
            %       - colony_size: the inital number of plants in the
            %         colony, which will then be increased up to colony_size as plants
            %         reproduce
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
            %        RES = studga(@myfunc,...)
            %     In this case, F = myfunc(X) returns the scalar function value F of
            %     the MYFUNC function evaluated at X.
            %
            %     FUN can also be an anonymous function:
            %        RES = studga(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
            %     returns X = [0;0].
            %
            %     clear classes,FUN=inline('sum(x.^2)'); n=100;
            %     lb=-20*ones(n,1); ub=-lb; x0=lb+(ub-lb).*rand(n,1);
            %     optimpot=struct('colony_size',20,'max_iter',10000,'max_time',180,'max_fcount',inf);
            %     RES=studga(@(x) FUN(x),x0,lb,ub,optimpot)
            
            % Reference: Dan Simon
            
            %   Copyright 2011 Junior Maih (junior.maih@gmail.com).
            %   $Revision: 7 $  $Date: 2011/05/27 17:23 $
            
            if nargin==0
                return
            end
            if nargin<5
                options=[];
            end
            % set these default options and change them if later on options is non-empty
            
            if ~isempty(options)
                all_prop=fieldnames(studga);
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
                    obj.f0=nan(1,n0);
                    for ii=1:n0
                        obj.f0(ii)=obj.Objective(obj.x0(:,ii),obj.vargs{:});
                    end
                    obj.fcount=obj.fcount+n0;
                else
                    obj.f0=obj.f0(1:n0);
                end
                [obj.f0,obj.x0]=resort(obj.f0,obj.x0);
                obj.best_fit=obj.f0(1);
                obj.best=obj.x0(:,1);
            end
            
            % Now find the peak
            obj.optimize;
        end
    end
end
