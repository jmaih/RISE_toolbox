classdef gampc < handle
    properties
        start_time
        stopping_created=false;
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
        best_fit
        xx
        ff
        number_of_parameters
        finish_time
        iter=0;
        fcount=0;
        optimizer='gampc';
    end
    methods(Static)
        function winners=tournament_selection(colony_size,choice_set)
            pool_size=3*colony_size; % pool size
            winners=nan(1,pool_size);
            for ii=1:pool_size
                % select the number of individuals to compete
                TcSize =randi(choice_set);
                % select the competitors
                randnum=nan(1,TcSize);
                for  tc=1:TcSize
                    randnum(tc) = randi(colony_size);
                end
                % It is assumed the population is sorted. In that case the
                % winner is the guy with the smallest index
                winners(ii) = min(randnum);
            end
        end
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
            %=============================
            arch_size=round(.5*obj.colony_size); % archive size
            %=============================
            [obj.xx(:,n0+1:end),obj.ff(n0+1:end),funevals]=...
                generate_candidates(obj.Objective,obj.lb,obj.ub,missing,...
                obj.restrictions,obj.penalty,obj.vargs{:});
            [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
            obj.fcount=obj.fcount+funevals;
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
                obj.iter=obj.iter+1;
                
                % Select the best performers
                tournament_bests=obj.tournament_selection(obj.colony_size,[2,3]);
                
                % Mutation
                offsprings=obj.mutation(tournament_bests);
                
                % Crossover
                offsprings=crossover(obj,offsprings,arch_size);
                
                obj.selection(offsprings,arch_size);
                
                % clear duplicates without re-evaluating the Objective function
                obj.xx=clear_duplicates(obj.xx,obj.lb,obj.ub,true);
                
                obj.memorize_best_solution;
                if rem(obj.iter,obj.verbose)==0 || obj.iter==1
                    restart=1;
                    fmin_iter=obj.best_fit;
                    disperse=dispersion(obj.xx(:,1:arch_size),obj.lb,obj.ub);
                    display_progress(restart,obj.iter,obj.best_fit,fmin_iter,...
                        disperse,obj.fcount,obj.optimizer);
                end
                stopflag=check_convergence(obj);
            end
            obj.finish_time=clock;
        end
        
        function offsprings=mutation(obj,tournament_bests)
            offsp_size=1:3:obj.colony_size;
            offsp_size=offsp_size(end)+2;
            offsprings=nan(obj.number_of_parameters,offsp_size);
            for ii=1:3:obj.colony_size
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
                    consec=unique([consec,randi(obj.colony_size)]);
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
            obj.fcount=obj.fcount+offsp_size;
            
            %  From both the archive individuals and the new offsprings, select
            %  the tournament_bests colony_size individuals for the next
            %  generation
            
            % 1-  Sort All polulation according to Objective value
            [all_fit,all_individuals] = resort(all_fit,all_individuals);
            
            % 2- copy the tournament_bests colony_size individuals into xx to start the new
            % generation
            obj.xx=all_individuals(:,1:obj.colony_size);
            obj.ff=all_fit(1:obj.colony_size);
            
%             [obj.xx,obj.ff,funevals]=rebuild_population(obj.xx,obj.ff,obj.Objective,obj.lb,obj.ub,0.03,obj.vargs{:});
%             [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
%             obj.fcount=obj.fcount+funevals;
        end
        function offsprings=crossover(obj,offsprings,arch_size)
            % Create an archive pool = 0.5*colony_size
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
            if isempty(obj.best_fit)||obj.ff(1)<obj.best_fit
                obj.best_fit=obj.ff(1);
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
            %   parameter vector), 'best_fit'(best function value), 'xx'(parameter
            %   vectors in the colony),'ff'(function values at xx)
            %   'max_iter','max_time','colony_size
            %
            %   RES = GAMPC(FUN,X0,[],LB,UB) starts at X0 and finds a minimum X to the
            %   function FUN, subject to the bound constraints LB and UB. FUN accepts
            %   input X and returns a scalar function value F evaluated at X. X0 may be
            %   a scalar or a vector.
            %
            %   RES = GAMPC(FUN,X0,[],LB,UB,OPTIONS) optimizes the function FUN under the
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
            %     optimpot=struct('colony_size',20,'max_iter',1000,'max_time',60,...
            %     'max_fcount',inf);
            %     RES=gampc(@(x) FUN(x),x0,[],lb,ub,optimpot)
            
            % Reference: Saber M. Elsayed, Ruhul A. Sarker and Daryl L.
            % Essam :" GA with a New Multi-Parent Crossover for Solving
            % IEEE-CEC2011 Competition Problems"
            
            %   Copyright 2011 Junior Maih (junior.maih@gmail.com).
            %   $Revision: 8 $  $Date: 2011/08/13 11:23 $
            
            if nargin==0
                return
            end
            if nargin<5
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
                obj.best_fit=obj.f0(1);
                obj.best=obj.x0(:,1);
            end
            
            % Now find the peak
            obj.optimize;
        end
    end
end