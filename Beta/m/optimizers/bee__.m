function [xbest,fbest,exitflag,output]=bee(objective,x0,lb,ub,options,varargin)

% this function sorts using f while the class alternative sorts using
% fitness. can this account for differences in results? in any case, the
% two functions seem to have fundamentally different behaviors. What about
% the generation of new bees?

%   Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   $Revision: 8 $  $Date: 2011/08/13 11:23 $

default_options = struct('MaxIter',1000,'MaxNodes',20,'MaxTime',3600,...
    'MaxFunEvals',inf,...
    'nonlcon',[],'rand_seed',100*sum(clock),...
    'verbose',10,'penalty',1e+8,...
    'max_trials',100,'start_time',clock,'stopping_created',false);

if nargin==0
    xbest=default_options;
    return
end
if nargin<5
    options=[];
end
ff=fieldnames(default_options);
for ifield=1:numel(ff)
    v=ff{ifield};
    if isfield(options,v)
        default_options.(v)=options.(v);
    end
end
optimizer='abc';
nonlcon=default_options.nonlcon;
if isempty(default_options.start_time)
    default_options.start_time=clock;
else
    if numel(default_options.start_time)~=6
        error([mfilename,':: wrong entry for start_time (should be same format as clock)'])
    end
end
verbose=default_options.verbose;
max_trials=default_options.max_trials;
MaxNodes=default_options.MaxNodes;
stopping_created=default_options.stopping_created;
rand_seed=default_options.rand_seed;
s = RandStream.create('mt19937ar','seed',rand_seed);
try
    RandStream.setGlobalStream(s);
catch %#ok<CTCH>
    RandStream.setDefaultStream(s); %#ok<SETRS>
end
output={'MaxIter','MaxTime','MaxFunEvals','start_time'};
not_needed=fieldnames(rmfield(default_options,output));
output=rmfield(default_options,not_needed);
output.funcCount = 0;
output.iterations = 0;
output.algorithm = optimizer;

%         [objective,xbest,options] = separateOptimStruct(objective);
npar = numel(x0);
food_number=round(.5*MaxNodes);

objective=fcnchk(objective,length(varargin));
pop = new_bees(x0);
for ib=2:MaxNodes
    pop(1,ib)=new_bees();
end
pop=sort_population(pop,output.iterations);
best_bee=[];
memorize_best_source();

% Now find the peak
optimize_bees();

% output
xbest=best_bee.x;
fbest=best_bee.f;
exitflag=1;

    function optimize_bees()
        % set and record the seed before we start drawing anything
        
        if ~stopping_created
            manual_stopping();
        end
        stopflag=check_convergence(output);
        while isempty(stopflag)
            output.iterations=output.iterations+1;
            send_employed_bees();
            send_onlooker_bees();
            memorize_best_source();
            send_scout_bees();
            if rem(output.iterations,verbose)==0 || output.iterations==1
                restart=1;
                fmin_iter=best_bee.f;
                disperse=dispersion([pop.x],lb,ub);
                display_progress(restart,output.iterations,best_bee.f,fmin_iter,...
                    disperse,output.funcCount,optimizer);
            end
            stopflag=check_convergence(output);
        end
        output.finish_time=clock;
    end

    function send_employed_bees()
        for ii=1:food_number
            pop(ii)=generate_mutant(ii);
        end
    end

    function send_onlooker_bees()
        fitness=[pop.fitness];
%         test=true;
%         if test
            probs= fitness/sum(fitness);
            cumprobs=[0,cumsum(probs)];
            cumprobs(end)=1;
            for t=1:food_number
                ii=find(cumprobs>rand,1,'first')-1;
                pop(ii)=generate_mutant(ii);
            end
%         else
%             prob=(0.9.*fitness./max(fitness))+0.1;
%             t=0;
%             ii=1;
%             while t<food_number
%                 if rand<prob(ii)
%                     t=t+1;
%                     pop(ii)=generate_mutant(ii);
%                 end
%                 ii=ii+1;
%                 if ii==food_number+1
%                     ii=1;
%                 end
%             end
%         end
    end

    function send_scout_bees()
        renew=find([pop.trial]>=max_trials);
        for ii=1:numel(renew)
            pop(renew(ii))=new_bees();
        end
    end

    function this=new_bees(x)
        if nargin==0
            x=lb+(ub-lb).*rand(npar,1);
        end
        this=wrapper(x);
        this.trial=0;
    end

    function memorize_best_source()
        sorted_bees=sort_population([pop,best_bee],output.iterations);
        best_bee=sorted_bees(1);
    end

    function mutant=generate_mutant(index)
        % generate a solution for index ii
        % a randomly chosen solution different from ii is used
        % for producing a mutant
        order=randperm(food_number); 
        order(order==index)=[];
        mutant=pop(index).x;
        
        bee_mutation=false;
        if bee_mutation
            donor_id=order(1);
            %     pick the parameters to change in the new solution
            change=min(fix(rand*npar)+1,npar);
            mutant(change)=mutant(change)+(mutant(change)-...
                pop(donor_id).x(change))*2.*(rand(numel(change),1)-.5);
        else
            if rand<.1,F1 = 0.5+0.3*randn;else F1 = 0.7+0.1*randn;end
            if rand<.1,F2 = 0.5+0.3*randn;else F2 = 0.7+0.1*randn;end
            F2=0;
%             F1=2*(rand-.5);
%             F2=2*(rand-.5);
            V=mutant+F1*(pop(order(2)).x-pop(order(3)).x)+F2*(pop(order(4)).x-pop(order(5)).x);
            % crossover
            Crmin=0.2;
            Crmax=0.7;
            Cr=Crmin+(Crmax-Crmin)*rand;
            crossers=rand(npar,1)<Cr;
            mutant(crossers)=V(crossers);
            %

        end
        mutant=new_bees(mutant);
        mutant=penalty_selection(mutant,pop(index),output.iterations);
        %         mutant=deb_selection(obj,ii,viol_strength,mutant,fit_mut,f_mut);
    end

    function this=wrapper(bird)
        if nargin==0
            this=evaluate_individual();
        else
            this=evaluate_individual(bird,objective,lb,ub,nonlcon,varargin{:});
            output.funcCount=output.funcCount+1;
        end
    end
end

function c=penalty_selection(a,b,k)
use_fitness=true;
c=a;
pa=dynamic_penalty(a.viol,k);
pb=dynamic_penalty(b.viol,k);
if use_fitness
    fpa=a.fitness-pa;
    fpb=b.fitness-pb;
    if fpb>fpa
        c=b;
    end
else
    ta=pa+a.f;
    tb=pb+b.f;
    if tb<ta
        c=b;
    end
end
end

