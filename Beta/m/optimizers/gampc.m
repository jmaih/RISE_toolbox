classdef gampc < handle
    properties
        start_time
        stopping_created=false;
        MaxNodes=20;
        lb
        ub
        x0
        f0
        MaxIter=1000;
        MaxTime=3600;
        MaxFunEvals=inf;
        rand_seed=100*sum(clock);
        penalty=1e+8;
        verbose=10;
        % algorithm specific options
        crossover_probability=0.01;
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
        best_fval
        xx
        ff
        number_of_parameters
        finish_time
        iterations=0;
        funcCount=0;
        optimizer='gampc';
    end
    methods(Static)
        function winners=tournament_selection(MaxNodes,choice_set)
            pool_size=3*MaxNodes; % pool size
            winners=nan(1,pool_size);
            for ii=1:pool_size
                % select the number of individuals to compete
                TcSize =randi(choice_set);
                % select the competitors
                randnum=nan(1,TcSize);
                for  tc=1:TcSize
                    randnum(tc) = randi(MaxNodes);
                end
                % It is assumed the population is sorted. In that case the
                % winner is the guy with the smallest index
                winners(ii) = min(randnum);
            end
        end
    end
    methods(Access=private)
        function obj=optimize(obj)
            obj.xx=nan(obj.number_of_parameters,obj.MaxNodes);
            obj.ff=nan(1,obj.MaxNodes);
            n0=size(obj.x0,2);
            if n0
                obj.ff(1:n0)=obj.f0;
                obj.xx(:,1:n0)=obj.x0(:,1:n0);
            end
            missing=obj.MaxNodes-n0;
            % set and record the seed before we start drawing anything
            s = RandStream.create('mt19937ar','seed',obj.rand_seed);
            RandStream.setDefaultStream(s);
            %=============================
            arch_size=round(.5*obj.MaxNodes); % archive size
            %=============================
            [obj.xx(:,n0+1:end),obj.ff(n0+1:end),~,funevals]=...
                generate_candidates(obj.Objective,obj.lb,obj.ub,missing,...
                obj.restrictions,obj.penalty,obj.vargs{:});
            [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
            obj.funcCount=obj.funcCount+funevals;
            obj.memorize_best_solution;
            
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
                
                % Select the best performers
                tournament_bests=obj.tournament_selection(obj.MaxNodes,[2,3]);
                
                % Mutation
                offsprings=obj.mutation(tournament_bests);
                
                % Crossover
                offsprings=crossover(obj,offsprings,arch_size);
                
                obj.selection(offsprings,arch_size);
                
                % clear duplicates without re-evaluating the Objective function
                obj.xx=clear_duplicates(obj.xx,obj.lb,obj.ub,true);
                
                obj.memorize_best_solution;
                if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
                    restart=1;
                    fmin_iter=obj.best_fval;
                    disperse=dispersion(obj.xx(:,1:arch_size),obj.lb,obj.ub);
                    display_progress(restart,obj.iterations,obj.best_fval,fmin_iter,...
                        disperse,obj.funcCount,obj.optimizer);
                end
                stopflag=check_convergence(obj);
            end
            obj.finish_time=clock;
        end
        
        function offsprings=mutation(obj,tournament_bests)
            offsp_size=1:3:obj.MaxNodes;
            offsp_size=offsp_size(end)+2;
            offsprings=nan(obj.number_of_parameters,offsp_size);
            for ii=1:3:obj.MaxNodes
                if rand<.1%
                    betta = 0.5+0.3*randn;
                else
                    betta = 0.7+0.1*randn;
                end
                
                consec=tournament_bests(ii:ii+2);
                % sort the selected three parents in ascending order, which is
                % equivalent to sorting them in according to their fitness
                consec= sort(consec);
                
                %%% Check the similarity between all selected individuals
                while ~isequal(numel(unique(consec)),3)
                    consec=unique([consec,randi(obj.MaxNodes)]);
                end
                offsprings(:,ii)=obj.xx(:,consec(1))+betta.*(obj.xx(:,consec(2))-obj.xx(:,consec(3)));
                offsprings(:,ii+1)=obj.xx(:,consec(2))+betta.*(obj.xx(:,consec(3))-obj.xx(:,consec(1)));
                offsprings(:,ii+2)=obj.xx(:,consec(3))+betta.*(obj.xx(:,consec(1))-obj.xx(:,consec(2)));
            end
            offsprings=recenter(offsprings,obj.lb,obj.ub);
        end
        function obj=selection(obj,offsprings,arch_size)
            % Group both elite and xx
            offsp_size=size(offsprings,2);
            all_individuals=[obj.xx(:,1:arch_size),offsprings];
            all_fit=[obj.ff(1:arch_size),nan(1,offsp_size)];
            
            % Calculated the fitness values for the neww offspring
            for ii=arch_size+(1:offsp_size)
                all_fit(ii)=obj.Objective(all_individuals(:,ii),obj.vargs{:});
            end
            obj.funcCount=obj.funcCount+offsp_size;
            
            %  From both the archive individuals and the new offsprings, select
            %  the tournament_bests MaxNodes individuals for the next
            %  generation
            
            % 1-  Sort All polulation according to Objective value
            [all_fit,all_individuals] = resort(all_fit,all_individuals);
            
            % 2- copy the tournament_bests MaxNodes individuals into xx to start the new
            % generation
            obj.xx=all_individuals(:,1:obj.MaxNodes);
            obj.ff=all_fit(1:obj.MaxNodes);
            
%             [obj.xx,obj.ff,funevals]=rebuild_population(obj.xx,obj.ff,obj.Objective,obj.lb,obj.ub,0.03,obj.vargs{:});
%             [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
%             obj.funcCount=obj.funcCount+funevals;
        end
        function offsprings=crossover(obj,offsprings,arch_size)
            % Create an archive pool = 0.5*MaxNodes
            archive=obj.xx(:,1:arch_size);
            % Randomized Operator
            for ii=1:size(offsprings,2)
                for jj=1:obj.number_of_parameters
                    if obj.crossover_probability>rand
                        pos= randi(arch_size);% select an individual from the archive pool
                        offsprings(jj,ii)=  archive(jj,pos); % exchange parameters
                    end
                end
            end
        end
        
        function obj=memorize_best_solution(obj)
            [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
            if isempty(obj.best_fval)||obj.ff(1)<obj.best_fval
                obj.best_fval=obj.ff(1);
                obj.best=obj.xx(:,1);
            end
        end
    end
    methods
        function obj=gampc(Objective,x0,f0,lb,ub,options,varargin)
            %  GAMPC attempts to find the global minimum of a constrained function of
            %  several variables.
            %   GAMPC attempts to solve problems of the form:
            %    min F(X)  subject to:  LB <= X <= UB   (bounds)
            %     X
            %
            %   RES = GAMPC constructs an object with the default optimization parameters.
            %   RES has fields that are intuitive to understand. 'lb' (lower bound),
            %   'ub' (upper bound),'x0',(vector of initial values),'f0'(function value
            %   at x0), 'vargs'(additional arguments of FUN),'penalty' (threshold
            %   function value beyond which the parameter draws are
            %   discarded),'Objective' (name of the Objective function), 'best'(best
            %   parameter vector), 'best_fval'(best function value), 'xx'(parameter
            %   vectors in the colony),'ff'(function values at xx)
            %   'MaxIter','MaxTime','MaxNodes
            %
            %   RES = GAMPC(FUN,X0,[],LB,UB) starts at X0 and finds a minimum X to the
            %   function FUN, subject to the bound constraints LB and UB. FUN accepts
            %   input X and returns a scalar function value F evaluated at X. X0 may be
            %   a scalar or a vector.
            %
            %   RES = GAMPC(FUN,X0,[],LB,UB,OPTIONS) optimizes the function FUN under the
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
            %        RES = gampc(@myfunc,...)
            %     In this case, F = myfunc(X) returns the scalar function value F of
            %     the MYFUNC function evaluated at X.
            %
            %     FUN can also be an anonymous function:
            %        RES = gampc(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[],[0 0])
            %     returns X = [0;0].
            %
            %     FUN=inline('sum(x.^2)'); n=100;
            %     lb=-20*ones(n,1); ub=-lb; x0=lb+(ub-lb).*rand(n,1);
            %     optimpot=struct('MaxNodes',20,'MaxIter',1000,'MaxTime',60,...
            %     'MaxFunEvals',inf);
            %     RES=gampc(@(x) FUN(x),x0,[],lb,ub,optimpot)
            
            % Reference: Saber M. Elsayed, Ruhul A. Sarker and Daryl L.
            % Essam :" GA with a New Multi-Parent Crossover for Solving
            % IEEE-CEC2011 Competition Problems"
            
            %   Copyright 2011 Junior Maih (junior.maih@gmail.com).
            %   $Revision: 8 $  $Date: 2011/08/13 11:23 $
            
            if nargin==0
                return
            end
            if nargin<6
                options=[];
            end
            
            if ~isempty(options)
                all_prop=fieldnames(gampc);
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
                if isempty(obj.f0)
                    obj.f0=nan(1,1:n0);
                    for ii=1:n0
                        obj.f0(ii)=obj.Objective(obj.x0(:,ii),obj.vargs{:});
                    end
                    obj.funcCount=obj.funcCount+n0;
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