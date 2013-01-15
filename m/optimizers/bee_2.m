function [xx,ff]=bee_2(Objective,x0,f0,lb,ub,options,varargin)

% Reference: Inspired from Karaboga

%   Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   $Revision: 8 $  $Date: 2011/08/13 11:23 $
if nargin<5
    options=[];
end

colony_size=[];
max_iter=[]; max_time=[]; max_fcount=[]; penalty=[];
rand_seed=[]; verbose=[];stopping_created=[];
known_optimum=[]; tol_fun=[];start_time=[];restrictions=[];

% optimizer-specific properties
max_trials=[];
fields={'colony_size','MaxNodes',50
    'max_iter','MaxIter',10000 % effectively max_iter/number_of_cycles
    'max_time','MaxTime',3600
    'max_fcount','MaxFunEvals',inf
    'rand_seed','rand_seed',100*sum(clock);
    'penalty','penalty',1e+8;
    'verbose','verbose',10
    'known_optimum','known_optimum',nan
    'tol_fun','TolFun',1e-6
    'start_time','start_time',clock
    'stopping_created','stopping_created',false
    'max_trials','max_trials',100
	'restrictions','restrictions',[]
    };
% change the items in the defaults
if ~isempty(options)
    options_prop=fieldnames(options);
    for ii=1:numel(options_prop)
        loc=find(strcmp(options_prop{ii},fields(:,1)));
        if isempty(loc)
            loc=find(strcmp(options_prop{ii},fields(:,2)));
        end
        if isempty(loc)
            error([mfilename,':: ',options_prop{ii},' is not a valid property'])
        end
        fields{loc,3}=options.(options_prop{ii});
    end
end

for ii=1:size(fields,1)
    eval([fields{ii,1},'=fields{ii,3};'])
end

if isempty(start_time)
    start_time=clock;
else
    if numel(start_time)~=6
        error([mfilename,':: wrong entry for start_time (should be same format as clock)'])
    end
end
number_of_parameters=size(lb,1);
Objective=fcnchk(Objective,length(varargin));
fixed=lb==ub;
variable=~fixed;
fixed_vals=lb(fixed);
npar_v=sum(variable);

% trim everything
lb=lb(variable);
ub=ub(variable);

obj=struct('known_optimum_reached',false,'max_fcount',max_fcount,...
    'fcount',0,'max_iter',max_iter,'iter',0,...
    'max_time',max_time,'start_time',start_time,'optimizer',mfilename);
n0=size(x0,2);

best_fval=[];
bestx=[];

if n0
    n0=min(n0,colony_size);
    x0=x0(variable,1:n0);
    if isempty(f0)
        f0=nan(1,1:n0);
        for ii=1:n0
            f0(ii)=evaluate(x0(:,ii));
        end
    else
        f0=f0(1:n0);
    end
    [f0,x0]=resort(f0,x0);
    best_fval=f0(1);
    bestx=x0(:,1);
end

optimizer=mfilename;
food_number=round(.5*colony_size);
xx=nan(npar_v,colony_size);
ff=nan(1,colony_size);
fitness=nan(1,colony_size);
trial=nan(1,colony_size);
n0=size(x0,2);
if n0
    fitness(1:n0)=compute_fitness(f0(1:n0));
    trial(1:n0)=0;
    ff(1:n0)=f0;
    xx(:,1:n0)=x0(:,1:n0);
end
missing=colony_size-n0;
% set and record the seed before we start drawing anything
s = RandStream.create('mt19937ar','seed',rand_seed);
RandStream.setDefaultStream(s);

% Now find the peak

[xx(:,n0+1:end),ff(n0+1:end),fitness(n0+1:end),...
    trial(n0+1:end)]=new_bees(missing);

memorize_best_source();

if ~stopping_created %#ok<*BDSCI,*BDLGI>
    manual_stopping();
else
    if ~exist('ManualStoppingFile.txt','file')
        manual_stopping();
    end
end
stopflag=check_convergence(obj);
while isempty(stopflag)
    obj.iter=obj.iter+1;
    send_employed_bees();
    send_onlooker_bees();
    memorize_best_source();
    proposal=[]; % proposal=xx(:,round(1+rand*(food_number-1)));
    send_scout_bees(proposal);
    if rem(obj.iter,verbose)==0 || obj.iter==1|| obj.iter==max_iter
        restart=1;
        fmin_iter=min(ff);
        disperse=dispersion(xx,lb,ub);
        display_progress(restart,obj.iter,best_fval,fmin_iter,...
            disperse,obj.fcount,optimizer);
    end
    if ~isnan(known_optimum) && abs(best_fval-known_optimum)<tol_fun
        obj.known_optimum_reached=true;
    end
    stopflag=check_convergence(obj);
end
xx=restore(xx);
bestx=restore(bestx);
xx(:,1)=bestx;
ff(1)=best_fval;

    function fy=evaluate(y)
        y=restore(y);
        fy=Objective(y,varargin{:});
        obj.fcount=obj.fcount+1;
    end

    function xx=restore(x)
        ncols=size(x,2);
        xx=nan(number_of_parameters,ncols);
        xx(fixed,:)=fixed_vals(:,ones(ncols,1));
        xx(variable,:)=x;
    end


    function send_scout_bees(xbest)
        renew=find(trial>=max_trials);
%         if ismember(1,renew)
%             [xstore,fstore]=local_reoptimize(@evaluate,xx(:,1),ff(:,1),lb,ub,[]);
%             if fstore(:,1)<ff(:,1)
%                 xx(:,1)=xstore(:,1);
%                 ff(:,1)=fstore(:,1);
%                 trial(1)=0;
%                 fitness(1)=compute_fitness(ff(:,1));
%                 renew=setdiff(renew,1);
%             end
%         end
        if nargin==1 && ~isempty(xbest)
            [xx(:,renew),ff(renew),fitness(renew),...
                trial(renew)]=new_bees4_renewal(xbest,numel(renew));
        else
            [xx(:,renew),ff(renew),fitness(renew),...
                trial(renew)]=new_bees(numel(renew));
        end
        [ff,xx,fitness,trial]=resort(ff,xx,fitness,trial);
    end
    function send_employed_bees()
        for i1=1:food_number
            generate_mutant(i1);
        end
    end
    function send_onlooker_bees()
        prob=(0.9.*fitness./max(fitness))+0.1;
        t=0;
        i1=1;
        while t<food_number
            if rand<prob(i1)
                t=t+1;
                generate_mutant(i1);
            end
            i1=i1+1;
            if i1==food_number+1
                i1=1;
            end
        end
    end
    function memorize_best_source
        [ff,xx,fitness,trial]=resort(ff,xx,fitness,trial);
        if isempty(best_fval)||ff(1)<best_fval
            best_fval=ff(1);
            bestx=xx(:,1);
        end
    end

    function [x,f,fit,trial]=new_bees(n)
        [x,f]=generate_candidates(@evaluate,lb,ub,n,restrictions,penalty);
        trial=zeros(1,n);
        fit=compute_fitness(f);
    end

    function [x,f,fit,trial]=new_bees4_renewal(xbest,n)
        x=bsxfun(@plus,lb,bsxfun(@times,ub-lb,rand(npar_v,n)));
        x=x+rand(1,1).*(xbest(:,ones(1,n))-x);
        x=recenter(x,lb,ub);
        f=nan(1,n);
        for i1=1:n
            f(i1)=evaluate(x(:,i1));
        end
        trial=zeros(1,n);
        fit=compute_fitness(f);
    end

    function generate_mutant(i1)
        % generate a solution for index i1
        % a randomly chosen solution different from i1 is used
        % for producing a mutant
        while 1
            donor_id=min(food_number,fix(rand*food_number)+1);
            if donor_id~=i1
                break
            end
        end
        mutant=xx(:,i1);
        
        change=round(1+rand*(npar_v-1));
        mutant(change)=mutant(change)+(mutant(change)-...
            xx(change,donor_id))*2.*(rand(numel(change),1)-.5);
        
%         mutant=mutant+(xx(:,donor_id)-mutant)*(1+randn);
        
        mutant=recenter(mutant,lb,ub);
        f_mut=evaluate(mutant);
        fit_mut=compute_fitness(f_mut);
        % we apply a greedy selection between i1 and the
        % mutant
        if fit_mut>fitness(i1)
            fitness(i1)=fit_mut;
            xx(:,i1)=mutant;
            ff(:,i1)=f_mut;
            trial(i1)=0;
        else
            trial(i1)=trial(i1)+1;
        end
    end

end
