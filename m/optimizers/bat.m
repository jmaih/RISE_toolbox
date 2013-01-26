function [xbest,fbest,exitflag,output]=bat(objective,x0,lb,ub,options,varargin)

% Reference: 
% Gaige Wang, Lihong Guo, Hong Duan, Luo Liu and Heqi Wang (2012): " A Bat
% Algorithm with Mutation for UCAV Path Planning". The Scientific World
% Journal, vol. 2012

default_options=struct('MaxIter',1000,'MaxFunEvals',150000,'TolX',1e-4,...
    'TolFun',1e-4,'nonlcon',[],'MaxNodes',20,'loudness',0.5,'pulse_rate',.5,...
    'min_frequency',0,'max_frequency',2,'alpha',.9,'gamma',.9);

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
MaxIter=default_options.MaxIter;  % Number of generations
MaxFunEvals=default_options.MaxFunEvals;
A=default_options.loudness;      % Loudness  (constant or decreasing)
r=default_options.pulse_rate;      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=default_options.min_frequency;         % Frequency minimum
Qmax=default_options.max_frequency;         % Frequency maximum
nonlcon=default_options.nonlcon;
alpha=default_options.alpha;
gamma=default_options.gamma;

output=struct('funcCount',0,'iterations',0);
output.funcCount=0;       % Total number of function evaluations
% Dimension of the search variables
d=numel(lb);           % Number of dimensions

% Initialize the population/solutions
bat_pop=wrapper(x0);
for i=2:n
    x=lb+(ub-lb).*rand(d,1);
    bat_pop(1,i)=wrapper(x);
end

bat_pop=sort_population(bat_pop,0);
% Find the initial best solution
best=bat_pop(1);

while output.iterations<MaxIter && output.funcCount<MaxFunEvals
    output.iterations=output.iterations+1;
    % Loop over all bats/solutions
    for i=1:n
        Qi=Qmin+(Qmin-Qmax)*rand;
        bat_pop(i).v=bat_pop(i).v+(bat_pop(i).x-best.x)*Qi;
        popi=wrapper(bat_pop(i).x+bat_pop(i).v);
        % Pulse rate
        if rand>bat_pop(i).r
%             r4=max(1,ceil(n*rand)); % select one of the best here. this should change
%             % The factor 0.001 limits the step sizes of random walks
%             Su=bat_pop(r4).x+0.001*randn(d,1);
            Su=best.x+0.001*randn(d,1);
            popu=wrapper(Su);
            popi=selection_process(popu,popi);
        end
        % generate a new solution by flying randomly
        epsilon=-1+(1--1)*rand;
        pop_rand=wrapper(bat_pop(i).x+epsilon*mean([bat_pop.A],2));
        popi=selection_process(pop_rand,popi);

        % Update if the solution improves, or not too loud
        choice=compare_individuals(popi,bat_pop(i));
        if choice==1 && rand<bat_pop(i).A
%         if choice==1 && rand<bat_pop(i).A
            bat_pop(i).x=popi.x;
            bat_pop(i).f=popi.f;
%             t=output.iterations;
%             bat_pop(i).r=r*(1-exp(-gamma*t));
%             bat_pop(i).A=alpha*bat_pop(i).A;
        end
        % Update the current best solution
        best=selection_process(best,popi);
    end
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', output.iterations, output.funcCount, best.f, mfilename)
end

xbest=best.x;
fbest=best.f;
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


