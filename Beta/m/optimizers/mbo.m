function [xbest,fbest,exitflag,output]=mbo(objective,x0,lb,ub,options,varargin)

% migrating birds optimization (MBO)

default_options = optimization_universal_options();

specific_options=struct('operator','de',...
    'iwo_sig_max',10/100,... % max standard deviation in IWO operator
    'iwo_sig_min',1/100,... % min standard deviation in IWO operator
    'm_',1,...% m_ : # of tours or wing flaps [10]
    'x_',1,...% x_ : # of neighbors solutions to be shared with the next solution [1]
    'k_',11);% k_ : # of neighbors solutions to be considered : speed of the flight [3]
%     'mutation_probability',0.01,...
%     'crossover_probability',1,...
%     'crossover_type',3,...

default_options=mergestructures(default_options,specific_options);
if nargin==0
    xbest=default_options;
    return
end

if nargin<5
    options=[];
end

[default_options]=mysetfield(default_options,options);

n_ =default_options.MaxNodes;
if rem(n_,2)==0
    n_=n_+1;
end
npar=numel(x0);
nonlcon=default_options.nonlcon;
m_=default_options.m_;
x_=default_options.x_;
k_=default_options.k_;
verbose=default_options.verbose;
operator=default_options.operator;
iwo_sig_max=default_options.iwo_sig_max;
iwo_sig_min=default_options.iwo_sig_min;

if iwo_sig_min>iwo_sig_max
    error('iwo_sig_min greater than iwo_sig_max')
end

if strcmp(operator,'iwo')
    search_range=ub-lb;
    % correct for inf
    search_range(isinf(search_range))=5;
    if any(search_range>10)
        warning('search range for some parameters exceeds 10...')
    end
    iwo_sig_min=iwo_sig_min*search_range;
    iwo_sig_max=iwo_sig_max*search_range;
end

output={'MaxIter','MaxTime','MaxFunEvals','start_time'};
not_needed=fieldnames(rmfield(default_options,output));
output=rmfield(default_options,not_needed);
output.funcCount = 0;
output.iterations = 0;
output.algorithm = [mfilename,'(',operator,')'];

[leader,Sleft,Sright]=initialize();
best_bird=leader;
% this could be small in principle and it is not defined in the paper/flowchart
p_=numel(Sleft);

Sleft_t_neighbor_set=wrapper();
Sright_t_neighbor_set=wrapper();

left_side=true;

if ~default_options.stopping_created
    manual_stopping();
end
stopflag=check_convergence(output);
while isempty(stopflag)
    output.iterations=output.iterations+1;
    if strcmp(operator,'iwo')
        coef=min((output.MaxIter-output.iterations)/output.MaxIter,(output.MaxFunEvals-output.funcCount)/output.MaxFunEvals);
        iwo_sig=iwo_sig_min+(iwo_sig_max-iwo_sig_min)*coef;
    end
    jj=0;
    while jj<m_
        jj=jj+1;
        improve_leader_solution();
        improve_other_solutions();
    end
    replace_leader();
    if rem(output.iterations,verbose)==0 || output.iterations==1
        restart=1;
        fmin_iter=best_bird.f;
        disperse=dispersion([pop.x],lb,ub);
        display_progress(restart,output.iterations,best_bird.f,fmin_iter,...
            disperse,output.funcCount,output.algorithm);
    end
    stopflag=check_convergence(output);
end

xbest=best_bird.x;
fbest=best_bird.f;
exitflag=1;

    function [leader,Sleft,Sright]=initialize()
        % randomly generate n_ initial solutions
        pop=wrapper(x0);
        for idraw=2:n_
            xdraw=lb+(ub-lb).*rand(npar,1);
            pop(idraw)=wrapper(xdraw);
        end
        % choose one of them as the leader
        viols=[pop.viol];
        fs=[pop.f];
        lead=find(viols==min(viols));
        lead=lead(fs(lead)==min(fs(lead)));
        leader=pop(lead);
        % build two lists: left and right
        % add remaining n_-1 solutions to the lists
        pop(lead)=[];
        Sleft=pop(1:2:n_-1);
        Sright=pop(2:2:n_-1);
    end
    function improve_leader_solution()
        % generate k_ neighbors of the leader
        nabo=generate_neighbors(leader,k_);
        % % %         ii=ii+k_;
        % sort the neighbors according to their fitness in non-decreasing
        nabo=sort_population(nabo);
        % order
        leader=selection_process(leader,nabo(1));
        
        nabo(1)=[];
        % add nb2, nb4, ... to neighbor set of Sleft(1)
        Sleft_t_neighbor_set=nabo(1:2:end);
        %         Sleft_t_neighbor_set=sort_population([Sleft_t_neighbor_set,nabo(1:2:end)]);
        % add nb3, nb5, ... to neighbor set of Sright(1)
        Sright_t_neighbor_set=nabo(2:2:end);
        %         Sright_t_neighbor_set=sort_population([Sright_t_neighbor_set,nabo(2:2:end)]);
        % which ones go out then?
        % update the best vector if possible
        best_bird=selection_process(best_bird,leader);
    end
    function improve_other_solutions()
        t=1;
        while t<p_
            % generate k_-x_ neighbors of slt and k_-x_ neighbors of srt
            nabo_left=generate_neighbors(Sleft(t),k_-x_,Sleft_t_neighbor_set);
            nabo_right=generate_neighbors(Sright(t),k_-x_,Sright_t_neighbor_set);
            % combine with x_ solutions from the front
            nabo_left=[nabo_left,Sleft_t_neighbor_set(1:x_)]; %#ok<*AGROW>
            nabo_right=[nabo_right,Sright_t_neighbor_set(1:x_)];
            % % %             ii=ii+2*(k_-x_);
            % sort the neighbors according to their fitness in a
            % non-decreasing manner
            nabo_left=sort_population(nabo_left);
            nabo_right=sort_population(nabo_right);
            Sleft(t)=selection_process(Sleft(t),nabo_left(1));
            Sright(t)=selection_process(Sright(t),nabo_right(1));
            % update the best vector if possible
            best_bird=selection_process(best_bird,Sleft(t));
            best_bird=selection_process(best_bird,Sright(t));
            % add nbl2, nbl3, ... to the neighborhood of sl(t+1)
            Sleft_t_neighbor_set=[Sleft_t_neighbor_set,nabo_left(2:end)];
            % add nbr2, nbr3, ... to the neighborhood of sr(t+1)
            Sright_t_neighbor_set=[Sright_t_neighbor_set,nabo_right(2:end)];
            t=t+1;
        end
    end
    function replace_leader()
        if left_side
            % move the leader to the end of the left list and asiwo_sign sl1 as
            % the new leader
            Sleft=[Sleft,leader];
            leader=Sleft(1);
            Sleft(1)=[];
        else
            Sright=[Sright,leader];
            leader=Sright(1);
            Sright(1)=[];
            % move the leader to the end of the right list and asiwo_sign sr1 as
            % the new leader
        end
        left_side=~left_side;
    end


    function mutants=generate_neighbors(stud,howmany,support_group)
        if nargin<3
            support_group=[];
        end
        nDonors=numel(support_group);
        % %     pick the parameters to change in the new solution
        mutants=wrapper();
        for dd=1:howmany
            if nDonors
                switch operator
                    case 'bee'
                        this=bee_operator(stud);
                    case 'de'
                        this=de_operator(stud);
                    case 'bbo'
                        this=bbo_operator(stud);
                    case 'iwo'
                        this=iwo_operator(stud);
                    case 'normal'
                        this=normal_operator(stud);
                    otherwise
                end
            else
                % drawing randomly will also help maintain diversity
                % here... but how to draw randomly without being wasteful?
                ss=randn(npar,1);
                this=stud.x+ss;
            end
            % % % %  this=recenter(this,lb,ub); this is
            % not needed since the vector is automatically adjusted before being
            % evaluated and stored
            mutants(1,dd)=wrapper(this);
        end
        
        function this=normal_operator(stud)
            this=stud.x+randn(npar,1);
        end
        function this=bee_operator(stud)
            this=stud.x;
            change=min(fix(rand*npar)+1,npar);
            nc=numel(change);
            donor_id=min(nDonors,fix(rand*nDonors)+1);
            donor=support_group(donor_id).x;
            this(change)=this(change)+(this(change)-donor(change))*2.*(rand(nc,1)-.5);
        end
        function this=de_operator(stud)
            F=0.7+0.1*randn;%(1+rand)*.5;
            Cr=.5;
            code='best_2';
            order_=randperm(nDonors);
            x15=[support_group(order_(1:5)).x];
            % mutation
            switch code
                case 'rand_1'
                    V=x15(:,1)+F*(x15(:,2)-x15(:,3));
                case 'rand_2'
                    V=x15(:,1)+F*(x15(:,2)-x15(:,3))+F*(x15(:,4)-x15(:,5));
                case 'best_1'
                    V=stud.x+F*(x15(:,2)-x15(:,3));
                case 'best_2'
                    V=stud.x+F*(x15(:,2)-x15(:,3))+F*(x15(:,4)-x15(:,5));
                case 'rand_best_1'
                    V=x15(:,1)+F*(stud.x-x15(:,3))+F*(x15(:,4)-x15(:,5));
                otherwise
            end
            % crossover
            this=stud.x;
            mutating_params=rand(npar,1)<Cr;
            this(mutating_params)=V(mutating_params);
        end
        function this=bbo_operator(stud)
            prob_mutate=0.01;
            [~,order]=sort_population([support_group,stud]);
            n=numel(order);
            k=1:n;
            II(order)=1-k/n;
            EE(order)=k/n;
            II(end)=[]; % the stud is never chosen
            im_probs=II/sum(II);
            cum_probs=[0,cumsum(im_probs)];
            cum_probs(end)=1;
            this=stud.x;
            migrating_params=rand(npar,1)<EE(end);
            for ipar=1:npar
                if migrating_params(ipar)
                    towards=find(cum_probs>rand,1,'first')-1;
                    this(ipar)=support_group(towards).x(ipar);
                end
                % gaussian mutation
                if prob_mutate>rand
                    this(ipar)=this(ipar)+randn;
                end
            end
        end
        function this=iwo_operator(stud)
            iwo_sig_eta=rand*iwo_sig;
            this=stud.x+iwo_sig_eta.*randn(npar,1);
        end
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

