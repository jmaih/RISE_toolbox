function [xbest,fbest,exitflag,output]=batm(objective,x0,lb,ub,options,varargin)

% Reference: 
% Gaige Wang, Lihong Guo, Hong Duan, Luo Liu and Heqi Wang (2012): " A Bat
% Algorithm with Mutation for UCAV Path Planning". The Scientific World
% Journal, vol. 2012

default_options = optimization_universal_options();

specific_options=struct('loudness',0.5,'pulse_rate',.5,...
    'min_frequency',0,'max_frequency',2,'alpha',.9,'gamma',.9,...
    'F',0.5);

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
A=default_options.loudness;      % Loudness  (constant or decreasing)
r=default_options.pulse_rate;      % Pulse rate (constant or decreasing)
F=default_options.F;      % Mutation probability
verbose=default_options.verbose;
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=default_options.min_frequency;         % Frequency minimum
Qmax=default_options.max_frequency;         % Frequency maximum
nonlcon=default_options.nonlcon;
alpha=default_options.alpha;
gamma=default_options.gamma;

output={'MaxIter','MaxTime','MaxFunEvals','start_time'};
not_needed=fieldnames(rmfield(default_options,output));
output=rmfield(default_options,not_needed);
output.funcCount = 0;
output.iterations = 0;
output.algorithm = mfilename;
% Dimension of the search variables
d=numel(lb);           % Number of dimensions

% Initialize the population/solutions
bat_pop=wrapper(x0);
for i=2:n
    x=lb+(ub-lb).*rand(d,1);
    bat_pop(1,i)=wrapper(x);
end

bat_pop=sort_population(bat_pop,0);
% Find the initial best_bat solution
best_bat=bat_pop(1);

if ~default_options.stopping_created
    manual_stopping();
end
stopflag=check_convergence(output);
while isempty(stopflag)
    output.iterations=output.iterations+1;
    % Loop over all bats/solutions
    for i=1:n
        r4=max(1,ceil(n*rand));
        Qi=Qmin+(Qmin-Qmax)*rand;
        vi=bat_pop(i).v+(bat_pop(i).x-best_bat.x)*Qi;
        Si=bat_pop(i).x+vi;
        popi=wrapper(Si);
        popi.v=vi;
        % Pulse rate
        if rand>bat_pop(i).r
            % The factor 0.001 limits the step sizes of random walks
            Su=best_bat.x+0.001*randn(d,1);
        else
            rr=randperm(n);
            Su=bat_pop(rr(1)).x+F*(bat_pop(rr(2)).x-bat_pop(rr(3)).x);
        end
        popu=wrapper(Su);
        kids=sort_population([popi,popu,bat_pop(r4)]);
        popk=kids(1);
        if rand<A && popk.f<bat_pop(r4).f
            bat_pop(r4)=popk;
        end
        % Update if the solution improves, or not too loud
        choice=compare_individuals(popi,bat_pop(i));
        if choice==1 || rand<bat_pop(i).A
%         if choice==1 && rand<bat_pop(i).A
            t=output.iterations;
            bat_pop(i)=popi;
            bat_pop(i).r=r*(1-exp(-gamma*t));
            bat_pop(i).A=alpha*bat_pop(i).A;
        end
    end
    % Update the current best_bat solution
    best_bat=selection_process(best_bat,popi);
    if rem(output.iterations,verbose)==0 || output.iterations==1
        restart=1;
        fmin_iter=best_bat.f;
        disperse=dispersion([bat_pop.x],lb,ub);
        display_progress(restart,output.iterations,best_bat.f,fmin_iter,...
            disperse,output.funcCount,output.algorithm);
    end
    stopflag=check_convergence(output);
end

xbest=best_bat.x;
fbest=best_bat.f;
exitflag=1;

    function this=wrapper(x)
        if nargin==0
            this=evaluate_individual();
            this.v=[];   % Velocities
            this.r=[];   % pulse rate
            this.A=[];
        else
            this=evaluate_individual(x,objective,lb,ub,nonlcon,varargin{:});
            this.v=0;   % Velocities
            this.r=r;   % pulse rate
            this.A =A;  % loudness
            output.funcCount=output.funcCount+1;
        end
    end

end


