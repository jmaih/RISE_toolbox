function [best,best_fval,obj]=gampc_2(Objective,x0,f0,lb,ub,Options,varargin)
%  GAMPC attempts to find the global minimum of a constrained function of
%  several variables.
%   GAMPC attempts to solve problems of the form:
%    min F(X)  subject to:  LB <= X <= UB   (bounds)
%     X

% Reference: Saber M. Elsayed, Ruhul A. Sarker and Daryl L.
% Essam :" GA with a New Multi-Parent Crossover for Solving
% IEEE-CEC2011 Competition Problems"

%   Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   $Revision: 8 $  $Date: 2011/08/13 11:23 $

colony_size=[];
max_iter=[];
% beta_0=[];
% beta_min=[];
% gamma=[];
start_time=[];
% alpha0=[];
max_time=[];
rand_seed=[];
penalty=[];
verbose=[];
max_fcount=[];
known_optimum=[];
crossover_probability=[];
stopping_created=[];
restrictions=[];

fields={'colony_size','MaxNodes',50
    'max_iter','MaxIter',10000 % effectively max_iter/number_of_cycles
    'max_time','MaxTime',3600
    'max_fcount','MaxFunEvals',inf
    'rand_seed','rand_seed',100*sum(clock);
    'penalty','penalty',1e+8;
    'verbose','verbose',10
    'start_time','start_time',clock
    'known_optimum','known_optimum',nan
    'crossover_probability','crossover_probability',.1
    'stopping_created','stopping_created',false
	'restrictions','restrictions',[]
    };

for ii=1:size(fields,1)
    if isfield(Options,fields{ii,1}) && ...
            ~isempty(Options.(fields{ii,1}))
        eval([fields{ii,1},'=Options.(fields{ii,1});'])
    elseif isfield(Options,fields{ii,2}) && ...
            ~isempty(Options.(fields{ii,2}))
        eval([fields{ii,1},'=Options.(fields{ii,2});'])
    else
        eval([fields{ii,1},'=fields{ii,3};'])
    end
end
if isempty(start_time)
    start_time=clock;
else
    if numel( start_time)~=6
        error([mfilename,':: wrong entry for start_time (should be same format as clock)'])
    end
end
number_of_parameters=size(lb,1);
Objective=fcnchk(Objective,length( varargin));
n0=size(x0,2);
obj.fcount=0;
if n0
    n0=min(n0, colony_size);
    x0= x0(:,1:n0);
    if isempty( f0)
        f0=nan(1,1:n0);
        for ii=1:n0
            f0(ii)= Objective( x0(:,ii), varargin{:});
        end
        obj.fcount= obj.fcount+n0;
    else
        f0= f0(1:n0);
    end
    [f0,x0]=resort( f0, x0);
end
xx=nan(number_of_parameters, colony_size);
ff=nan(1, colony_size);
n0=size( x0,2);
if n0
    ff(1:n0)= f0;
    xx(:,1:n0)= x0(:,1:n0);
end
missing= colony_size-n0;
% set and record the seed before we start drawing anything
s = RandStream.create('mt19937ar','seed', rand_seed);
RandStream.setDefaultStream(s);
%=============================
arch_size=round(.5*colony_size); % archive size
%=============================
[ xx(:,n0+1:end), ff(n0+1:end),funevals]=...
    generate_candidates( Objective, lb, ub,missing,...
    restrictions,penalty, varargin{:});
[ff,xx]=resort(ff,xx);
best=xx(:,1);
best_fval=ff(1);
obj.fcount= obj.fcount+funevals;
pool_size=3*colony_size; % pool size
offsp_size=1:3: colony_size;
offsp_size=offsp_size(end)+2;

memorize_best_solution;

if ~ stopping_created %#ok<BDSCI,BDLGI>
    manual_stopping;
else
    if ~exist('ManualStoppingFile.txt','file')
        manual_stopping;
    end
end
obj.iter=0;
obj.max_iter=max_iter;
obj.start_time=start_time;
obj.max_time=max_time;
obj.max_fcount=max_fcount;
obj.optimizer=mfilename;
stopflag=check_convergence(obj);
while isempty(stopflag)
    obj.iter= obj.iter+1;
    
    % Select the best performers
    tournament_bests= tournament_selection([2,3]);
    
    % Mutation
    offsprings= mutation(tournament_bests);
    
    % Crossover
    offsprings=crossover(offsprings);
    
    selection(offsprings);
    
    memorize_best_solution;
    if rem( obj.iter, verbose)==0 ||  obj.iter==1
        restart=1;
        fmin_iter= best_fval;
        disperse=dispersion( xx(:,1:arch_size),lb,ub);
        display_progress(restart, obj.iter, best_fval,fmin_iter,...
            disperse, obj.fcount,mfilename);
    end
    if ~isnan(known_optimum) && abs(best_fval-known_optimum)<1e-8
        obj.known_optimum_reached=true;
    end
    stopflag=check_convergence(obj);
end
obj.finish_time=clock;

    function winners=tournament_selection(choice_set)
        winners=nan(1,pool_size);
        for i1=1:pool_size
            % select the number of individuals to compete
            TcSize =randi(choice_set);
            % select the competitors
            randnum=nan(1,TcSize);
            for  tc=1:TcSize
                randnum(tc) = randi(colony_size);
            end
            % It is assumed the population is sorted. In that case the
            % winner is the guy with the smallest index
            winners(i1) = min(randnum);
        end
    end

    function offsprings=mutation(tournament_bests)
        offsprings=nan( number_of_parameters,offsp_size);
        for i1=1:3: colony_size
            if rand<1%
                betta = 0.5+0.3*randn;
            else
                betta = 0.7+0.1*randn;
            end
            
            consec=tournament_bests(i1:i1+2);
            % sort the selected three parents in ascending order, which is
            % equivalent to sorting them in according to their fitness
            consec= sort(consec);
            
            %%% Check the similarity between all selected individuals
            while ~isequal(numel(unique(consec)),3)
                consec=unique([consec,randi(colony_size)]);
            end
            test=0;
            switch test
                case 0
                    offsprings(:,i1)=xx(:,consec(1))+betta.*(xx(:,consec(2))- xx(:,consec(3)));
                    offsprings(:,i1+1)=xx(:,consec(2))+betta.*(xx(:,consec(3))- xx(:,consec(1)));
                    offsprings(:,i1+2)=xx(:,consec(3))+betta.*(xx(:,consec(1))- xx(:,consec(2)));
                case 1
                    u1=rand(number_of_parameters,1);
                    o1=u1<.33; o2=u1>=.33&u1<.66; o3=~o1&~o2;
                    offsprings(o1,i1)=xx(o1,consec(1));
                    offsprings(o2,i1)=xx(o2,consec(2));
                    offsprings(o3,i1)=xx(o3,consec(3));
                    
                    offsprings(o1,i1+1)=xx(o1,consec(2));
                    offsprings(o2,i1+1)=xx(o2,consec(3));
                    offsprings(o3,i1+1)=xx(o3,consec(1));
                    
                    offsprings(o1,i1+2)=xx(o1,consec(3));
                    offsprings(o2,i1+2)=xx(o2,consec(1));
                    offsprings(o3,i1+2)=xx(o3,consec(2));
                case 2
                    offsprings(:,i1)=xx(:,consec(1));
                    id=round(1+rand*(number_of_parameters-1));
                    offsprings(id,i1)=offsprings(id,i1)+(2*rand-1)*(xx(id,consec(2))- xx(id,consec(3)));
                    
                    offsprings(:,i1+1)=xx(:,consec(2));
                    id=round(1+rand*(number_of_parameters-1));
                    offsprings(id,i1+1)=offsprings(id,i1+1)+(2*rand-1)*(xx(id,consec(3))- xx(id,consec(1)));
                    
                    offsprings(:,i1+2)=xx(:,consec(3));
                    id=round(1+rand*(number_of_parameters-1));
                    offsprings(id,i1+2)=offsprings(id,i1+2)+(2*rand-1)*(xx(id,consec(1))- xx(id,consec(2)));
                case 3
                    SD=.5*(max(xx(:,consec),[],2)-min(xx(:,consec),[],2));
                    offsprings(:,i1)=xx(:,consec(1))+SD.*randn(number_of_parameters,1);
                    offsprings(:,i1+1)=xx(:,consec(2))+SD.*randn(number_of_parameters,1);
                    offsprings(:,i1+2)=xx(:,consec(3))+SD.*randn(number_of_parameters,1);
                otherwise
            end
        end
        offsprings=recenter(offsprings,lb,ub);
    end

    function selection(offsprings)
        % Group both elite and xx
        all_individuals=[ xx(:,1:arch_size),offsprings];
        all_fit=[ff(1:arch_size),nan(1,offsp_size)];
        
        % Calculate the fitness values for the neww offspring, replacing
        % the duplicates
        [all_individuals,all_fit,funevals]=...
            rebuild_population(all_individuals,all_fit,Objective,...
            lb,ub,0.005,varargin{:});
        obj.fcount= obj.fcount+funevals;
        
        %  From both the archive individuals and the new offsprings, select
        %  the tournament_bests colony_size individuals for the next
        %  generation
        
        % 1-  Sort All according to Objective value
        [all_fit,all_individuals] = resort(all_fit,all_individuals);
        
        % 2- copy the tournament_bests colony_size individuals into xx to start the new
        % generation
        xx=all_individuals(:,1:colony_size);
        ff=all_fit(1:colony_size);
        [xx,ff]=rebuild_population(xx,ff,Objective,lb,ub,0.01,varargin{:});
        
%         [ff,xx]=resort(ff,xx);
%         obj.fcount= obj.fcount+funevals;

    end
    function offsprings=crossover(offsprings)
        % Create an archive pool = 0.5*colony_size
        archive= xx(:,1:arch_size);
        % Randomized Operator
        for i1=1:size(offsprings,2)
            for jj=1: number_of_parameters
                if  crossover_probability>rand
                    pos= randi(arch_size);% select an individual from the archive pool
                    offsprings(jj,i1)=  archive(jj,pos); % exchange parameters
                end
            end
        end
    end

    function memorize_best_solution()
        [ff,xx]=resort(ff,xx);
        if isempty(best_fval)|| ff(1)<best_fval
            best_fval= ff(1);
            best= xx(:,1);
        end
    end
end