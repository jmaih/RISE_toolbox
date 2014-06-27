function [xbest,fbest,exitflag,output]=diff_lev_bat(objective,x0,lb,ub,options,varargin)

% Reference:
% Jian Xie, Yongquan Zhou, Huan Chen (201?): A Novel Bat Algorithm Based on
% Differential Operator and Levy-Flights Trajectory

default_options = optimization_universal_options();

specific_options=struct('loudness_0_min',1,'loudness_0_max',2,...
    'pulse_rate_min',0,'pulse_rate_max',0.1,'Fmin',0,'Fmax',1,...
    'alpha',.9);

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
A0_min=default_options.loudness_0_min;      % Loudness  (constant or decreasing)
A0_max=default_options.loudness_0_max;      % Loudness  (constant or decreasing)
rmin=default_options.pulse_rate_min;      % Pulse rate (constant or decreasing)
rmax=default_options.pulse_rate_max;      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
nonlcon=default_options.nonlcon;
alpha=default_options.alpha;
verbose=default_options.verbose;
Fmin=default_options.Fmin;
Fmax=default_options.Fmax;
lambdamin=1+sqrt(eps);
lambdamax=3;


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
%--------------------------------------------------------------------------
disperse=inf;
search_range=ub-lb;
search_range(isinf(search_range))=5;
sig=2/100*search_range;
%--------------------------------------------------------------------------

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
    t=output.iterations;
    coef=max(t/output.MaxIter,output.funcCount/output.MaxFunEvals);
    % Loop over all bats/solutions
    for b=1:n
        if disperse<1e-3
            x=best_bat.x+sig.*randn(npar,1);
%             x=lb+(ub-lb).*rand(npar,1);
            bat_b=wrapper(x);
            bat_pop(b).x=bat_b.x;
            bat_pop(b).f=bat_b.f;
        else
            b1=rand;
            b2=rand;
            mu=rand;
            lambda=lambdamin+(lambdamax-lambdamin)*rand;
            F1=((Fmin-Fmax)*coef+Fmax)*b1;
            F2=(Fmin+coef*(Fmax-Fmin))*b2;
            order=randperm(n);
            xr=[bat_pop(order(1:4)).x];
            xi=best_bat.x+F1*(xr(:,1)-xr(:,2))+F2*(xr(:,3)-xr(:,4));
            xi=xi+mu*sign(rand(npar,1)-.5)*t^(-lambda);
            popi=wrapper(xi);
            % Pulse rate (if high enough)
            if bat_pop(b).r>rand
                epsilon=2*rand(npar,1)-1;
                Su=best_bat.x+epsilon*mean([bat_pop.A],2);
                popu=wrapper(Su);
                popi=selection_process(popu,popi);
            end
            if bat_pop(b).A>rand
                eta=2*rand(npar,1)-1;
                pop_rand=wrapper(best_bat.x+eta*mean([bat_pop.r],2));
                popi=selection_process(pop_rand,popi);
            end
            
            choice=compare_individuals(popi,bat_pop(b));
            if choice==1
                bat_pop(b).x=popi.x;
                bat_pop(b).f=popi.f;
                choice=compare_individuals(bat_pop(b),best_bat);
                if choice==1
                    bat_pop(b).A=alpha*bat_pop(b).A;
                    bat_pop(b).r=bat_pop(b).r*coef^3;
                end
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
            this.r=[];   % pulse rate
            this.A=[];
        else
            this=evaluate_individual(x,objective,lb,ub,nonlcon,varargin{:});
            this.r=rmin+(rmax-rmin)*rand;   % pulse rate
            this.A =A0_min+(A0_max-A0_min)*rand;
            output.funcCount=output.funcCount+1;
        end
    end

end


