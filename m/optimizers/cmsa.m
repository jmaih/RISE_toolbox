classdef cmsa < handle
    properties
        TolFun=1e-6;
        known_optimum=nan;
        stopping_created=false;
        start_time
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
        mu_prop=.25; % must be between (0,1)
        sig0=5;
        sig_min=0.1;
        init_var=1;
        max_stagnate=200;
		restrictions
    end
    %     properties(Constant, Hidden = true)
    %     end
    properties (SetAccess = private, Hidden = true)
        known_optimum_reached=false;
        vargs={};
    end
    properties(SetAccess=protected)
        Objective
        best
        best_fval
        restart_best_fit
        xx
        ff
        number_of_parameters
        finish_time
        iterations=0;
        funcCount=0;
        optimizer='cmsa';
        Cov
    end
    methods(Static)
        function [mm,C,sig_rec]=updates(C,xx_search,ss,weights,sig,tau_c)
            SST=0;
            sig_rec=0;
            mm=0;
            for ii=1:numel(weights)
                mm=mm+weights(ii)*xx_search(:,ii); %
                SST=SST+weights(ii)*(ss(:,ii)*ss(:,ii)');
                sig_rec=sig_rec+weights(ii)*sig(ii);
            end
            C=(1-1/tau_c)*C+1/tau_c*SST;
            C=.5*(C+C');
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
            ss = RandStream.create('mt19937ar','seed',obj.rand_seed);
            RandStream.setDefaultStream(ss);
            [obj.xx(:,n0+1:end),obj.ff(n0+1:end),funevals]=...
                generate_candidates(obj.Objective,obj.lb,obj.ub,missing,...
                obj.restrictions,obj.penalty,obj.vargs{:});
            [obj.ff,obj.xx]=resort(obj.ff,obj.xx);
            obj.funcCount=obj.funcCount+funevals;
            obj.memorize_best_solution;
            
            %=========================
            mu=round(obj.mu_prop*obj.MaxNodes);
            ii=1:mu;
            weights=log(mu+.5)-log(ii);
            weights=weights/sum(weights);
            tau_c=1+obj.number_of_parameters*(obj.number_of_parameters+1)/(2*mu); % time constant
            sig_rec=obj.sig0;
            C=diag(obj.init_var*ones(obj.number_of_parameters,1));
            
            mm=0;
            for ii=1:mu
                mm=mm+weights(ii)*obj.xx(:,ii); %
            end
            stagnate=0;
            restart=1;
            obj.restart_best_fit=obj.best_fval;
            %=========================
            
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
                
                break_it=false;
                if stagnate==obj.max_stagnate
                    [C,mm,sig_rec,obj]=re_start(obj);
                    stagnate=0;
                    restart=restart+1;
                    if isempty(mm) % I could not find a stud
                        disp([mfilename,':: STOPPING UNDER UNUSUAL CIRCUMSTANCE OF FAILED RESTART'])
                        break_it=true;
                    end
                end
                if break_it
                    break
                end
                [ff_search,xx_search,ss,sig,obj.funcCount]=recombination(obj,C,mm,sig_rec);
                
                if ff_search(1)<obj.best_fval
                    % record the best population so far
                    obj.ff=ff_search;
                    obj.xx=xx_search;
                    obj.Cov=C;
                end
                if ff_search(1)<obj.restart_best_fit
                    obj.restart_best_fit=ff_search(1);
                    stagnate=0;
                else
                    stagnate=stagnate+1;
                end

                xx_search=clear_duplicates(xx_search,obj.lb,obj.ub);
                
                [mm,C,sig_rec]=obj.updates(C,xx_search,ss,weights,sig,tau_c);
                
                obj.memorize_best_solution;
                if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
                    disperse=dispersion(xx_search,obj.lb,obj.ub);
                    display_progress(restart,obj.iterations,obj.best_fval,ff_search(1),...
                        disperse,obj.funcCount,obj.optimizer);
                end
                if ~isnan(obj.known_optimum) && abs(obj.best_fval-obj.known_optimum)<obj.TolFun
                    obj.known_optimum_reached=true;
                end
                stopflag=check_convergence(obj);
            end
            obj.finish_time=clock;
        end
        
        function [C,mm,sig_rec,obj]=re_start(obj)
            % try generating a stud
            iter0=0;
            success=false;
            mm=[];
            while iter0<20 && ~success
                iter0=iter0+1;
                try %#ok<TRYNC>
                    [mm,fm,funevals]=generate_candidates(...
                        obj.Objective,obj.lb,obj.ub,obj.MaxNodes,...
						obj.restrictions,obj.penalty,obj.vargs{:});
                    obj.funcCount=obj.funcCount+funevals;
                    success=true;
                    [fm,mm]=resort(fm,mm);
                    mm=mm(:,1);
                    obj.restart_best_fit=fm(1);
                end
            end
            sig_rec=obj.sig0;
            C=diag(obj.init_var*ones(obj.number_of_parameters,1));
            %             C=obj.Cov;
        end
        
        function [ff_search,xx_search,ss,sig,funcCount]=recombination(obj,C,mm,sig_rec)
            % recombination for the cmsa guys
            tau=1/sqrt(2*obj.number_of_parameters); % learning parameter
            sig=max(obj.sig_min,sig_rec)*exp(tau*randn(1,obj.MaxNodes));
            CC=chol(C,'lower');
            ss=CC*randn(obj.number_of_parameters,obj.MaxNodes);
            z=bsxfun(@times,sig,ss);
            xx_search=mm(:,ones(1,obj.MaxNodes))+z;
            
            xx_search=recenter(xx_search,obj.lb,obj.ub);
            ff_search=xx_search(1,:);
            for l=1:obj.MaxNodes
                ff_search(l)=obj.Objective(xx_search(:,l),obj.vargs{:});
            end
            funcCount=obj.funcCount+obj.MaxNodes;
            [ff_search,xx_search,ss,sig]=resort(ff_search,xx_search,ss,sig);
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
        function obj=cmsa(Objective,x0,f0,lb,ub,options,varargin)
            %  CMSA attempts to find the global minimum of a constrained function of
            %  several variables.
            %   CMSA attempts to solve problems of the form:
            %    min F(X)  subject to:  LB <= X <= UB   (bounds)
            %     X
            %
            %   RES = CMSA constructs an object with the default optimization parameters.
            %   RES has fields that are intuitive to understand. 'lb' (lower bound),
            %   'ub' (upper bound),'x0',(vector of initial values),'f0'(function value
            %   at x0), 'vargs'(additional arguments of FUN),'penalty' (threshold
            %   function value beyond which the parameter draws are
            %   discarded),'Objective' (name of the Objective function), 'best'(best
            %   parameter vector), 'best_fval'(best function value), 'xx'(parameter
            %   vectors in the colony),'ff'(function values at xx)
            %   'MaxIter','MaxTime','MaxNodes
            %
            %   RES = CMSA(FUN,X0,[],LB,UB) starts at X0 and finds a minimum X to the
            %   function FUN, subject to the bound constraints LB and UB. FUN accepts
            %   input X and returns a scalar function value F evaluated at X. X0 may be
            %   a scalar or a vector.
            %
            %   RES = CMSA(FUN,X0,[],LB,UB,OPTIONS) optimizes the function FUN under the
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
            %        RES = cmsa(@myfunc,...)
            %     In this case, F = myfunc(X) returns the scalar function value F of
            %     the MYFUNC function evaluated at X.
            %
            %     FUN can also be an anonymous function:
            %        RES = cmsa(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[],[0 0])
            %     returns X = [0;0].
            %
            %     FUN=inline('sum(x.^2)'); n=100;
            %     lb=-20*ones(n,1); ub=-lb; x0=lb+(ub-lb).*rand(n,1);
            %     optimpot=struct('MaxNodes',20,'MaxIter',1000,'MaxTime',60,...
            %     'MaxFunEvals',inf);
            %     RES=cmsa(@(x) FUN(x),x0,[],lb,ub,optimpot)
            
            % Reference: Inspired from Beyer and Sendhoff (????):
            % "Covariance Matrix Adaptation Revisited -the CMSA Evolution
            % Strategy-
            
            %   Copyright 2011 Junior Maih (junior.maih@gmail.com).
            %   $Revision: 8 $  $Date: 2011/08/13 11:23 $
            
            if nargin==0
                return
            end
            if nargin<4
                options=[];
            end
            
            if ~isempty(options)
                all_prop=fieldnames(cmsa);
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