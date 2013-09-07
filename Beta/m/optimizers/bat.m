function [xbest,fbest,exitflag,output]=bat(objective,x0,lb,ub,options,varargin)

% Reference: 
% Gaige Wang, Lihong Guo, Hong Duan, Luo Liu and Heqi Wang (2012): " A Bat
% Algorithm with Mutation for UCAV Path Planning". The Scientific World
% Journal, vol. 2012

default_options = optimization_universal_options();

specific_options=struct('loudness_0',0.5,'loudness_min',0.5,'pulse_rate',.5,...
    'variable_loudness',true,'variable_pulse_rate',false,'min_frequency',0,'max_frequency',100,...
    'alpha',.9,'gamma',.9);

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
A0=default_options.loudness_0;      % Loudness  (constant or decreasing)
Amin=default_options.loudness_min;
r=default_options.pulse_rate;      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
freq_min=default_options.min_frequency;         % Frequency minimum
freq_max=default_options.max_frequency;         % Frequency maximum
nonlcon=default_options.nonlcon;
alpha=default_options.alpha;
gamma=default_options.gamma;
verbose=default_options.verbose;
variable_loudness=default_options.variable_loudness;
variable_pulse_rate=default_options.variable_pulse_rate;

output={'MaxIter','MaxTime','MaxFunEvals','start_time'};
not_needed=fieldnames(rmfield(default_options,output));
output=rmfield(default_options,not_needed);
output.funcCount = 0;
output.iterations = 0;
output.algorithm = mfilename;
% Dimension of the search variables
npar=numel(lb);           % Number of dimensions

% Initialize the population/solutions
bat_pop=wrapper(x0);
for b=2:n
    x=lb+(ub-lb).*rand(npar,1);
    bat_pop(1,b)=wrapper(x);
end

% % probabilities
% cprobs=n:-1:1;
% cprobs=[0,cumsum(cprobs/sum(cprobs))];
% cprobs(end)=1;

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
    for b=1:n
        freq_i=freq_min+(freq_min-freq_max)*rand;
        bat_pop(b).v=bat_pop(b).v+(bat_pop(b).x-best_bat.x)*freq_i;
        popi=wrapper(bat_pop(b).x+bat_pop(b).v);
        % Pulse rate (if high enough)
        if bat_pop(b).r>rand
            %             Su=pop(loc).x+0.001*randn(npar,1);
%             loc=find(cprobs>rand,1,'first')-1;
            epsilon=2*rand(npar,1)-1;
            Su=best_bat.x+epsilon*mean([bat_pop.A],2);
%             Su=bat_pop(loc).x+epsilon*mean([bat_pop.A],2);
            popu=wrapper(Su);
            popi=selection_process(popu,popi);
        end
        % generate a new solution by flying randomly
%         pop_rand=wrapper(best_bat.x+0.001*randn(npar,1));
        pop_rand=wrapper(bat_pop(b).x+0.001*randn(npar,1));
        popi=selection_process(pop_rand,popi);
        
        % Update if the solution improves, or not too loud
        choice=compare_individuals(popi,bat_pop(b));
        %         if choice==1 && rand<bat_pop(b).A
        if choice==1 || rand<bat_pop(b).A
            bat_pop(b).x=popi.x;
            bat_pop(b).f=popi.f;
            t=output.iterations;
            if variable_loudness
                bat_pop(b).A=alpha*bat_pop(b).A;
            end
            bat_pop(b).A=max(Amin,bat_pop(b).A);
            if variable_pulse_rate
                bat_pop(b).r=bat_pop(b).r*(1-exp(-gamma*t))/(1-exp(-gamma*(t-1)));
            end
        end
        % Update the current best_bat solution
        best_bat=selection_process(best_bat,popi);
    end
    if rem(output.iterations,verbose)==0 || output.iterations==1
        restart=1;
        fmin_iter=best_bat.f;
        disperse=dispersion([bat_pop.x],lb,ub);
        display_progress(restart,output.iterations,best_bat.f,fmin_iter,...
            disperse,output.funcCount,output.algorithm);
    end
    stopflag=check_convergence(output);
    bat_pop=sort_population(bat_pop,output.iterations);
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
            if variable_pulse_rate
                this.r=rand;   % pulse rate
            else
                this.r=r;   % pulse rate
            end
                this.A =A0;  % loudness_0
%             if variable_loudness
%                 this.A =Amin+(Amax-Amin)*rand;  % loudness_0
%             else
%             end
            output.funcCount=output.funcCount+1;
        end
    end

end


