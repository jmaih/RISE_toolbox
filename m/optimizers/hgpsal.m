function [xbest,fbest,exitflag,output]=hgpsal(objective,x0,lb,ub,options,varargin)

default_options = struct('MaxIter',300,'MaxFunEvals',150000,...
    'TolX',1e-4,'TolFun',1e-4,...
    'nonlcon',[],...
    'lambda_min',-1e+12,'lambda_max',1e+12,'delta_max',1e+12,...
    'mu_min',1e-12,'gamma_',0.5,'epsilon_star',1e-12,'eta_star',1e-6,...
    'mu',1,'eta',1,'pie',0.5,'tau',0.5);

% If just 'defaults' passed in, return the default options in X
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
% load the options
lambda_min=default_options.lambda_min;
lambda_max=default_options.lambda_max;
delta_max=default_options.delta_max;
mu_min=default_options.mu_min;
gamma_=default_options.gamma_;
epsilon_star=default_options.epsilon_star;
eta_star=default_options.eta_star;
mu=default_options.mu;
eta=default_options.eta;
pie=default_options.pie;
tau=default_options.tau;
MaxIter=default_options.MaxIter;
nonlcon=default_options.nonlcon;

output.funcCount=0;
output.iterations=1;
mm=0;
lambda=zeros(mm,1); % in min max, i=1,...,mm
cx=zeros(mm,1); % eguality constraints
delta=0; % in [0,delta_max], i=1:pp
this=wrapper(x0);
gx=this.viol; % inequality constraints
% pp=numel(this.viol);
npar=numel(lb);

% 1-compute epsilon
epsil=epsilon();
E=inf;

while output.iterations<MaxIter && E>eta_star && epsil>epsilon_star
    % find and ej-global minimizer xj of subproblem (2) using algorithm 1
    % so that
    hybrid_genetic_pattern_search();
    
    gx=this.viol;
    delta=max(0,min(delta+gx/mu,delta_max)); % i=1:pp
    
    E=big_e();
    if E<=eta
        lambda=lambda+max(lambda_min,min(cx/mu,lambda_max)); % i=1:mm
    else
        mu=max(mu_min,gamma_*mu);
    end
    eta=pie*eta;
    epsil=epsilon();
    output.iterations=output.iterations+1;
end

xbest=this.x;
fbest=this.f;
exitflag=1;


    function e=epsilon()
        e=max(epsilon_star,tau/(1+norm(lambda)+norm(delta)+1/mu));
    end
    function hybrid_genetic_pattern_search()
        this=genetic_algorithm(this);
        this=hooke_jeeves(this);
    end


    function obj=wrapper(bird)
        if nargin==0
            obj=evaluate_individual();
        else
            obj=evaluate_individual(bird,objective,lb,ub,nonlcon,varargin{:});
            % add the penalty
            obj.f=big_phi();
            output.funcCount=output.funcCount+1;
        end
        function phie=big_phi()
            gx=obj.viol; % inequality constraints
            def=delta+gx/mu;
            def=norm(def(def>0));
            phie=obj.f+lambda'*cx+1/(2*mu)*norm(cx)^2+0.5*mu*(def-norm(delta));
        end
    end

    function E=big_e()
        E=max([norm(cx),norm(gx),max(delta.*abs(gx))])/(1+norm(this.x));
    end

    function this=genetic_algorithm(this)
        pc=0.9;
        kmax=200;
        ss=20;
        ee=2;
        eta_c=20;
        pm=1/npar;
        eta_m=20;
        pop=wrapper(this.x);
        for ipop=2:ss
            pop(1,ipop)=wrapper(lb+(ub-lb).*rand(npar,1));
        end
        
        k=0;
        while k<kmax;
            pop=sort_population(pop,k);
            % select by tournamentss ss-ee points
            children=pop(ee+1:end);
            for ichild=1:0.5*(ss-ee)
                choice=tournament_selection(2);
                w=simulated_binary_crossover([pop(choice).x]);
                for aa=1:size(w,2)
                    w(:,aa)=polynomial_mutation(w(:,aa));
                    children((ichild-1)*2+aa)=wrapper(w(:,aa));
                end
            end
            pop=[pop,children]; %#ok<AGROW>
            % elitism
            pop(ee+1:end)=sort_population(pop(ee+1:end),k);
            pop=pop(1:ss);
            k=k+1;
            fprintf('iter %8.4g   fval  %8.4g   f-count  %8.12g\n',k,pop(1).f,output.funcCount);
        end
        keyboard
                
        function choice=tournament_selection(n)
            choice=[];
            order=randperm(ss);
            for ii=1:n
                indices=order(1:2);
                tmp=sort_population(pop(indices));
                if isequal(tmp(1),pop(indices(1)))
                    hit=indices(1);
                else
                    hit=indices(2);
                end
                choice=[choice,hit]; %#ok<AGROW>
                order(order==hit)=[];
                tmp=randperm(ss-ii);
                order=order(tmp);
            end
        end
        % select z(1) and z(2) randomly from the pool
                
        function t=polynomial_mutation(w)
            t=w;
            mutants=rand(npar,1)<pm;
            for ipar=1:npar
                if ~mutants(ipar)
                    continue
                end
                ri=rand;
                if ri<.5
                    iota=(2*ri)^(1/(1+eta_m))-1;
                else
                    iota=1-(2*(1-ri))^(1/(1+eta_m));
                end
                t(ipar)=w(ipar)+(ub(ipar)-lb(ipar))*iota;
            end
        end
        
        function w=simulated_binary_crossover(z)
            % with probability pc create two new points
            w=z;
            crossers=rand(npar,1)<pc;
            for ipar=1:npar
                if ~crossers(ipar)
                    continue
                end
                ri=rand;
                if ri<=.5
                    b=(2*ri)^(1/(1+eta_c));
                else
                    b=(1/(2*(1-ri)))^(1/(1+eta_c));
                end
                w(ipar,1)=0.5*((1+b)*z(ipar,1)+(1-b)*z(ipar,2));
                w(ipar,2)=0.5*((1-b)*z(ipar,1)+(1+b)*z(ipar,2));
            end
        end
    end
end