function [xbest,fbest,exitflag,output]=iwo(objective,x0,lb,ub,options,varargin)
    
% IWO: Invasive Weed Optimization
%   No Detailed explanation yet

% Reference: 
% Abhronil Sengupta, Tathagata Chakraborti, Amit Konar and Atulya K. Nagar
% (????): An Intelligent Invasive Weed Optimization: a Q-learning approach

default_options = optimization_universal_options();

specific_options=struct('Smin',0,'Smax',5,'sig_max',10/100,'sig_min',1/100,...
    'krigging_threshold',0);

default_options=mergestructures(default_options,specific_options);

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
n=default_options.MaxNodes;      % Population size, typically 10 to 40
nonlcon=default_options.nonlcon;
verbose=default_options.verbose;
Smin=default_options.Smin;
Smax=default_options.Smax;
sig_max=default_options.sig_max;
sig_min=default_options.sig_min;
krigging_threshold=default_options.krigging_threshold;

if sig_min>sig_max
    error('sig_min greater than sig_max')
end

search_range=ub-lb;
% correct for inf
search_range(isinf(search_range))=5;
if any(search_range>10)
    warning('search range for some parameters exceeds 10...')
end
sig_min=sig_min*search_range;
sig_max=sig_max*search_range;

output={'MaxIter','MaxTime','MaxFunEvals','start_time'};
not_needed=fieldnames(rmfield(default_options,output));
output=rmfield(default_options,not_needed);
output.funcCount = 0;
output.iterations = 0;
output.algorithm = mfilename;
% Dimension of the search variables
npar=numel(lb);           % Number of dimensions

% Initialize the population/solutions
pop=wrapper(x0);
for i=2:n
    x=lb+(ub-lb).*rand(npar,1);
    pop(1,i)=wrapper(x);
end

pop=sort_population(pop,output.iterations);
% Find the initial best_weed solution
best_weed=pop(1);

if ~default_options.stopping_created
    manual_stopping();
end
stopflag=check_convergence(output);
MaxIter=output.MaxIter;
MaxFunEvals=output.MaxFunEvals;
krigging_threshold=(1-krigging_threshold);
while isempty(stopflag)
    output.iterations=output.iterations+1;
    % Loop over all bats/solutions
    coef=min((MaxIter-output.iterations)/MaxIter,(MaxFunEvals-output.funcCount)/MaxFunEvals);
    sig=sig_min+(sig_max-sig_min)*coef;
    best_fit=max([pop.fitness]);
    worst_fit=min([pop.fitness]);
    denom=best_fit-worst_fit;
    seeds=wrapper();
    krigging=coef>=krigging_threshold;
    for i=1:n
        nseeds=ceil(Smin+(Smax-Smin)*(best_fit-pop(i).fitness)/denom);
        seeds_i=wrapper();
        for jseed=1:nseeds
            sig_eta=rand*sig;
            this=pop(i).x+sig_eta.*randn(npar,1);
            seeds_i(1,jseed)=wrapper(this);
        end
        if krigging
            seeds_i=sort_population([seeds_i,pop(i)]);
            pop(i)=seeds_i(1);
        else
            seeds=[seeds,seeds_i];
        end
    end
    pop=sort_population([pop,seeds],output.iterations);
    if numel(pop)>n
        pop=pop(1:n);
    end
    best_weed=pop(1);
    if rem(output.iterations,verbose)==0 || output.iterations==1
        restart=1;
        fmin_iter=best_weed.f;
        disperse=dispersion([pop.x],lb,ub);
        display_progress(restart,output.iterations,best_weed.f,fmin_iter,...
            disperse,output.funcCount,output.algorithm);
    end
    stopflag=check_convergence(output);
end

xbest=best_weed.x;
fbest=best_weed.f;
exitflag=1;

    function this=wrapper(x)
        if nargin==0
            this=evaluate_individual();
        else
            this=evaluate_individual(x,objective,lb,ub,nonlcon,varargin{:});
            output.funcCount=output.funcCount+1;
        end
    end

end



