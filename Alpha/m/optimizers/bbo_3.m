function [xx,ff,obj]=bbo_3(Objective,x0,f0,lb,ub,Options,varargin)

stopping_created=[];start_time=[];MaxNodes=[];
MaxIter=[];MaxTime=[];MaxFunEvals=[];rand_seed=[];penalty=[];
verbose=[]; TolFun=[]; migration_model=[]; maximum_immigration_rate=[];
maximum_emigration_rate=[]; mutation_probability=[];
habitat_modification_probability=[]; elitism=[]; mutation_type=[];
known_optimum=[]; restrictions=[];
Fields={'stopping_created',false
    'start_time',clock
    'MaxNodes',20
    'MaxIter',1000
    'MaxTime',3600
    'MaxFunEvals',inf
    'rand_seed',100*sum(clock)
    'penalty',1e+8
    'verbose',10
    'migration_model',3
    'maximum_immigration_rate',1
    'maximum_emigration_rate',1
    'mutation_probability',0.005
    'habitat_modification_probability',1
    'elitism',2
    'mutation_type','gaussian' % 'cauchy','uniform'
    'known_optimum',nan
    'TolFun',1e-6
	'restrictions',[]
    };
for ii=1:size(Fields,1)
    if isfield(Options,Fields{ii,1})
        eval([Fields{ii,1},'=Options.(Fields{ii,1});'])
    else
        eval([Fields{ii,1},'=Fields{ii,2};'])
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

obj=struct('known_optimum_reached',false,'MaxFunEvals',MaxFunEvals,...
    'funcCount',0,'MaxIter',MaxIter,'iterations',0,'finish_time',[],...
    'MaxTime',MaxTime,'start_time',start_time,'optimizer',mfilename);
n0=size(x0,2);
if n0
    n0=min(n0,MaxNodes);
    x0=x0(variable,1:n0);
    if isempty(f0)
        f0=nan(1,1:n0);
        for ii=1:n0
            f0(ii)=evaluate(x0(:,ii));
        end
    else
        f0=f0(1:n0);
    end
end
xx=nan(npar_v,MaxNodes);
ff=nan(1,MaxNodes);
if n0
    ff(1:n0)=f0;
    xx(:,1:n0)=x0(:,1:n0);
end
missing=MaxNodes-n0;
s = RandStream.create('mt19937ar','seed',rand_seed);
RandStream.setDefaultStream(s);
[xx(:,n0+1:end),ff(n0+1:end)]=...
    generate_candidates(@evaluate,lb,ub,missing,...
    restrictions,penalty);
[ff,xx]=resort(ff,xx);
best_fval=ff(1);

% Now find the peak
if ~stopping_created %#ok<*BDSCI,*BDLGI>
    manual_stopping();
else
    if ~exist('ManualStoppingFile.txt','file')
        manual_stopping();
    end
end

[immigration_rate,emigration_rate]=migration_rates();
memorize_best_solution();

stopflag=check_convergence(obj);
while isempty(stopflag)
    obj.iterations=obj.iterations+1;
    islands=migration();
    mutation();
    xx=clear_duplicates(xx,lb,ub); %clear_duplicates;
    selection();
    memorize_best_solution();
    if rem(obj.iterations,verbose)==0 || obj.iterations==1
        restart=1;
        fmin_iter=best_fval;
        disperse=dispersion(xx,lb,ub);
        display_progress(restart,obj.iterations,best_fval,fmin_iter,...
            disperse,obj.funcCount,mfilename);
    end
    if ~isnan(known_optimum) && abs(best_fval-known_optimum)<TolFun
        obj.known_optimum_reached=true;
    end
    stopflag=check_convergence(obj);
end
xx=restore(xx);
obj.finish_time=clock;

    function fy=evaluate(y)
        y=restore(y);
        fy=Objective(y,varargin{:});
        obj.funcCount=obj.funcCount+1;
    end

    function xx=restore(x)
        ncols=size(x,2);
        xx=nan(number_of_parameters,ncols);
        xx(fixed,:)=fixed_vals(:,ones(ncols,1));
        xx(variable,:)=x;
    end

    function memorize_best_solution()
        [ff,xx]=resort(ff,xx);
        if isempty(best_fval)||ff(1)<best_fval
            best_fval=ff(1);
            %             best=xx(:,1);
        end
    end


    function mutation()
        % Mutate only the worst half of the solutions
        kstart=round(MaxNodes/2);
        k=MaxNodes-kstart+1;
        tmp=islands(:,kstart:end);
        %             mutants=bsxfun(@lt,rand(number_of_parameters,k),mutation_rates(kstart:end));
        mutants=mutation_probability>rand(npar_v,k);
        switch mutation_type
            case 'uniform'
                mutations=bsxfun(@plus,lb,...
                    bsxfun(@times,rand(npar_v,k),ub-lb));
                mutations=mutations(mutants);
            case 'gaussian'
                number_of_mutations=sum(sum(mutants));
                mutations=tmp(mutants)+randn(number_of_mutations,1);
            case 'cauchy'
                Onesk=ones(1,k);
                mutations=truncated_cauchy_inv(rand(npar_v,k),0,1,lb(:,Onesk),ub(:,Onesk));
                mutations=mutations(mutants);
                %                 case 'levy'
            otherwise
                error([mfilename,':: mutation type ',mutation_type,' not implemented'])
        end
        tmp(mutants)=mutations;
        islands(:,kstart:end)=tmp;
        
        islands=recenter(islands,lb,ub);
    end

    function islands=migration()
        % Compute immigration rate and extinction rate for each species count.
        % lambda(i) is the immigration rate for individual i.
        % mu(i) is the extinction rate for individual i.
        % This routine assumes the population is sorted from most fit to least fit.
        mu=emigration_rate;
        lambda=immigration_rate;
        % select habitats to modify
        modify=habitat_modification_probability>rand(1,MaxNodes);
        lambda_scale = (lambda-min(lambda))/(max(lambda)-min(lambda));
        % initialize so that those that are not modified remain
        % untouched
        islands=xx;
        cmu=mu/sum(mu);
        cmu=cumsum(cmu);cmu(end)=1;
        for kk=1:MaxNodes
            if modify(kk)
                for jj = 1:npar_v
                    if rand<lambda_scale(kk)
                        % Pick a habitat from which to obtain a feature
                        SelectIndex=find(cmu>rand,1,'first');
                        islands(jj,kk) = xx(jj,SelectIndex);
                    end
                end
            end
        end
    end

    function [immigration_rate,emigration_rate]=migration_rates()
        I=maximum_immigration_rate;
        E=maximum_emigration_rate;
        n=MaxNodes;
        
        k=n:-1:1; % <--- k=1:n;
        switch migration_model
            case 1 % Constant immigration and linear emigration model
                Lambda=.5*I*ones(1,n);
                Mu=k/n*E;
            case 2 % Linear immigration and constant emigration model
                Lambda=I*(1-k/n);
                Mu=.5*E*ones(1,n);
            case 3 % Linear migration model
                Lambda=I*(1-(k-1)/n);
                Mu=(k-1)/n*E;
            case 4 % Trapezoidal migration model
                i1=ceil((n+1)/2);
                Lambda=nan(1,n);
                Mu=nan(1,n);
                Lambda(1:i1)=I;
                Mu(1:i1)=2*E/n*k(1:i1);
                Lambda(i1+1:end)=2*I*(1-k(i1+1:end)/n);
                Mu(i1+1:end)=E;
            case 5 % Quadratic migration model
                Lambda=I*(k/n-1).^2;
                Mu=E*(k/n).^2;
            case 6 % Sinusoidal migration model
                Lambda=I/2*(cos(pi*k/n)+1);
                Mu=E/2*(-cos(pi*k/n)+1);
        end
        immigration_rate=Lambda;
        emigration_rate=Mu;
    end

    function selection()
        HSI=ff;
        new_habitats=find(any(xx-islands));
        cc=numel(new_habitats);
        for i1=1:cc
            id=new_habitats(i1);
            HSI(id)=evaluate(islands(:,id));
        end
        [HSI,islands]=resort(HSI,islands);
        % bring in the elites
        index=0;
        for i1=1:elitism
            contains_elite=any(sum(bsxfun(@minus,islands,xx(:,i1)).^2,1))==0;
            if ~contains_elite
                index=index+1;
                islands(:,end-index+1)=xx(:,i1);
                HSI(end-index+1)=ff(i1);
            end
        end
        [ff,xx]=resort(HSI,islands);
        % this step is redone in memorize_best_solution below. But I do
        % it here just to be tidy. The philosophy is that the best guy
        % is always first.
    end

end