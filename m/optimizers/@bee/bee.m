classdef bee %< handle
    % BEE global optimization of constrained functions using bees
    %
    % Constructor
    % ------------
    %
    % - [bee](bee/bee)
    
    properties
        stopping_created=false;
        start_time
        MaxNodes=20
        lb
        ub
        x0
        f0
        violation_strength_0
        fitness_0
        iterations=0
        funcCount=0
        MaxIter=1000
        MaxTime=3600
        MaxFunEvals=inf
        rand_seed=[]
        penalty=1e+8
        verbose=10
        % optimizer-specific properties
        max_trials=100
        nonlcon
        max_genes_change = 1 % number of genes change
    end
    %     properties(Constant, Hidden = true)
    %     end
    properties (SetAccess = private, Hidden = true)
        food_number
        trial
        fitness
        vargs={};
        restrictions_options = struct('restrictions_in_objective',false,...
            'returns_retcode',false,'restrictions_same_weights',true,...
            'allow_restrictions_violations',true) % restrictions are treated by Deb
        bank % holds the candidates if they exceed what is required...
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
            %  BEE attempts to find the global minimum of a constrained function of several variables.
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
            
            nonlinear_constraint_in_objective=isfield(options,'nonlcon') &&...
                ~isempty(options.nonlcon) && islogical(options.nonlcon);
            
            if ~isfield(options,'nonlcon')||isempty(options.nonlcon)||...
                    nonlinear_constraint_in_objective
                
                options.nonlcon=@(z)[lb-z;z-ub];
                
            else
                
                options.nonlcon=@(z)[options.nonlcon(z,varargin{:});lb-z;z-ub];
                
            end
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
            
            obj.restrictions_options.restrictions_in_objective=nonlinear_constraint_in_objective;
            
            obj.lb=lb;
            
            obj.ub=ub;
            
            obj.x0=x0;
            
            obj.f0=real(f0);
            
            obj.violation_strength_0=imag(f0);
            
            if isempty(obj.start_time)
                
                obj.start_time=clock;
                
            else
                
                if numel(obj.start_time)~=6
                    
                    error([mfilename,':: wrong entry for start_time (should be same format as clock)'])
                    
                end
                
            end
            
            obj.vargs=varargin;
            
            obj.number_of_parameters=size(obj.lb,1);
            
            if ischar(Objective)
                
                Objective=str2func(Objective);
                
            end
            
            obj.Objective=@newObjective;
            
            n0=size(obj.x0,2);
            
            if isempty(obj.f0)
                
                obj.violation_strength_0=zeros(1,n0);
                
                obj.f0=nan(1,n0);
                
                for ii=1:n0
                    
                    [obj.f0(ii),obj.violation_strength_0(ii)]=...
                        utils.estim.eval_objective_and_restrictions(...
                        obj.x0(:,ii),obj.Objective,obj.nonlcon,...
                        obj.restrictions_options,obj.vargs{:});
                    
                end
                
                obj.funcCount=obj.funcCount+n0;
                
            end
            
            obj.fitness_0=utils.optim.compute_fitness(obj.f0);
            
            [obj.f0,obj.x0,obj.fitness_0,~,obj.violation_strength_0]=sort_bees_detail(...
                obj.f0,obj.x0,obj.fitness_0,zeros(1,n0),obj.violation_strength_0);
            
            obj=memorize_best_source(obj);
            
            m=obj.MaxNodes;
            
            obj.bank=struct('x',{},'f',{},'violation_strength',{});
            
            if n0>m
                
                for k=1:n0-m
                    
                    obj.bank(k).x=obj.x0(:,k+m);
                    
                    obj.bank(k).f=obj.f0(k+m);
                    
                    obj.bank(k).violation_strength=obj.violation_strength_0(k+m);
                    
                end
                
                obj.x0=obj.x0(:,1:m);
                
                obj.f0=obj.f0(1:m);
                
                obj.violation_strength_0=obj.violation_strength_0(1:m);
                
                obj.fitness_0=obj.fitness_0(1:m);
                
            end
            
            % Now find the peak
            obj=optimize_bees(obj);
            
            function [f,retcode,viol]=newObjective(varargin)
                
                retcode=0;
                
                if obj.restrictions_options.restrictions_in_objective
                    
                    [f,retcode,viol]=Objective(varargin{:});
                    
                else
                    
                    f=Objective(varargin{:});
                    
                    viol=[];
                    
                end
                
                viol=[viol;obj.nonlcon(varargin{1})];
                
            end
            
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
    
    obj.violation_strength(1:n0)=obj.violation_strength_0(1:n0);
    
end

% save memory
obj.fitness_0=obj.fitness_0(1);
obj.f0=obj.f0(1);
obj.x0=obj.x0(:,1);
obj.violation_strength_0=obj.violation_strength_0(1);

missing=obj.MaxNodes-n0;
% set and record the seed before we start drawing anything

if ~isempty(obj.rand_seed)
    
    rng(obj.rand_seed)
    
end

[obj.xx(:,n0+1:end),obj.ff(n0+1:end),...
    obj.violation_strength(n0+1:end),funevals,...
    obj.fitness(n0+1:end),obj.trial(n0+1:end)]=...
    new_bees(obj.Objective,obj.lb,obj.ub,missing,obj.nonlcon,...
    obj.restrictions_options,obj.penalty,obj.max_trials,obj.vargs{:});

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
    
    obj=send_employed_bees(obj);
    
    obj=send_onlooker_bees(obj);
    
    obj=memorize_best_source(obj);
    
    obj=send_scout_bees(obj);
    
    utils.optim.display_progress(obj);
    
    stopflag=utils.optim.check_convergence(obj);
    
end

obj.finish_time=clock;

end

function obj=send_scout_bees(obj)

renew=find(obj.trial>=obj.max_trials);

% first look up the bank
n=numel(obj.bank);

while n>0 && ~isempty(renew)
   
    pinga=randi(n);
    
    pingb=renew(1);
    
    obj.xx(:,pingb)=obj.bank(pinga).x;
    
    obj.ff(pingb)=obj.bank(pinga).f;
    
    obj.violation_strength(pingb)=obj.bank(pinga).violation_strength;
    
    obj.fitness(pingb)=utils.optim.compute_fitness(obj.ff(pingb));
    
    obj.trial(pingb)=0;
    
    renew=renew(2:end);
    
    obj.bank(pinga)=[];
    
    n=n-1;
    
end

if isempty(renew)
    
    return
    
end

[obj.xx(:,renew),obj.ff(renew),obj.violation_strength(renew),funevals,...
    obj.fitness(renew),obj.trial(renew)]=new_bees(obj.Objective,obj.lb,...
    obj.ub,numel(renew),obj.nonlcon,obj.restrictions_options,...
    obj.penalty,obj.max_trials,obj.vargs{:});

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
    
    tmp=deb_selection(obj,1,obj.best_viol_strength,obj.best,...
        obj.best_fitness,obj.best_fval);
    
    tmpviol=tmp.violation_strength(1);
    
    tmpf=tmp.ff(1);
    
    tmpx=tmp.xx(:,1);
    
    if isempty(obj.best_fval)||...
            (tmpviol<obj.best_viol_strength) || ...
            (tmpviol==obj.best_viol_strength && tmpf<obj.best_fval)
        
        obj.best_fval=tmpf;
        
        obj.best=tmpx;
        
        obj.best_viol_strength=tmpviol;
        
    end
    
end

end

function obj=sort_bees(obj)

[obj.ff,obj.xx,obj.fitness,obj.trial,obj.violation_strength]=sort_bees_detail(...
    obj.ff,obj.xx,obj.fitness,obj.trial,obj.violation_strength);

end

function [FF,XX,FIT,TRIALS,VIOLS]=sort_bees_detail(FF0,XX0,FIT0,TRIALS0,VIOLS0)
VS=unique(VIOLS0);

FF=nan(size(FF0));

XX=nan(size(XX0));

FIT=nan(size(FIT0));

TRIALS=nan(size(TRIALS0));

VIOLS=nan(size(VIOLS0));

iter0=0;

for iloc=1:numel(VS)
    
    this=VIOLS0==VS(iloc);
    
    nthis=sum(this); % this is logical and so the sum has to be taken
    
    [FF(iter0+(1:nthis)),XX(:,iter0+(1:nthis)),FIT(iter0+(1:nthis)),...
        TRIALS(iter0+(1:nthis)),VIOLS(iter0+(1:nthis))]=...
        utils.optim.resort(FF0(this),XX0(:,this),FIT0(this),...
        TRIALS0(this),VIOLS0(this));
    
    iter0=iter0+nthis;
    
end

end

function [x,f,viol_strength,funevals,fit,trial]=new_bees(objective,lb,ub,n,...
    nonlcon,opt,penalty,max_trials,varargin)

[x,f,viol_strength,funevals]=utils.optim.generate_candidates(objective,lb,...
    ub,n,max_trials,nonlcon,opt,penalty,varargin{:});

if f>=penalty
    
    viol_strength = penalty;
    
end

trial=zeros(1,n);

fit=utils.optim.compute_fitness(f);

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
nch=min(randi(ceil(0.5*obj.number_of_parameters)),obj.max_genes_change);

change=min(fix(rand(nch,1)*obj.number_of_parameters)+1,obj.number_of_parameters);

mutant(change)=mutant(change)+(mutant(change)-...
    obj.xx(change,donor_id))*2.*(rand(nch,1)-.5);

mutant(change)=utils.optim.recenter(mutant(change),obj.lb(change),obj.ub(change));

[f_mut,viol_strength]=utils.estim.eval_objective_and_restrictions(mutant,...
    obj.Objective,obj.nonlcon,obj.restrictions_options);

if f_mut>=obj.penalty
    
    viol_strength=obj.penalty;
    
end

obj.funcCount=obj.funcCount+1;

fit_mut=utils.optim.compute_fitness(f_mut);

obj=deb_selection(obj,ii,viol_strength,mutant,fit_mut,f_mut);

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