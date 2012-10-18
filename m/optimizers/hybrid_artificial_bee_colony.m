function [xbest,fbest,obj]=hybrid_artificial_bee_colony(objective,x0,lb,ub,Options)

if nargin<5
    Options=struct('max_time',3600/60*5);
end
stopping_created=[]; start_time=[]; colony_size=[]; max_iter=[]; 
max_time=[]; max_fcount=[]; verbose=[]; abc=[];limit=[]; 
F=[]; cr=[]; debug=[]; threshold=[];
Fields={'stopping_created',false
    'start_time',clock
    'colony_size',20
    'max_iter',10000
    'max_time',3600
    'max_fcount',inf
    'verbose',10
    'abc',3
    'limit',100
    'F',0.5 % F=2*rand; %[0,2]
    'cr',0.8 % F=2*rand; %[0,2]
    'debug',false
    'threshold',1e-5;
    };
abc=ceil(abc);
if abc<1
    error([mfilename,':: abc=[{1}|2|3]'])
    % 1: basic abc
    % 2: abc with DE-like employed bee step
    % 3: abc with DE-like employed bee step + DE step
end
de_step=abc>2;
for ii=1:size(Fields,1)
    if isfield(Options,Fields{ii,1})
        eval([Fields{ii,1},'=Options.(Fields{ii,1});'])
    else
        eval([Fields{ii,1},'=Fields{ii,2};'])
    end
end
n_empl=colony_size;

obj=struct('known_optimum_reached',false,'max_fcount',max_fcount,...
    'fcount',0,'max_iter',max_iter,'iter',0,'finish_time',[],...
    'max_time',max_time,'start_time',start_time,'optimizer',mfilename);

n_onlook=ceil(0.5*n_empl);
nparticles=ceil(0.75*n_empl);
trail=zeros(1,n_empl);
npar=numel(lb);
missing=n_empl-size(x0,2);
x=[x0,bsxfun(@plus,lb,bsxfun(@times,ub-lb,rand(npar,missing)))];
f=nan(1,n_empl);
% % % % fitness=nan(1,n_empl);
for ii=1:n_empl
    f(ii)=evaluate(x(:,ii));
end
xbest=[];
fbest=inf;
memorize_best_food_source();

% Now find the peak
if ~stopping_created %#ok<*BDSCI,*BDLGI>
    manual_stopping();
else
    if ~exist('ManualStoppingFile.txt','file')
        manual_stopping();
    end
end
stopflag=check_convergence(obj);

while  isempty(stopflag)
    obj.iter=obj.iter+1;
    
    % -1) search new food source and evaluate its quality
    for ie=1:n_empl
        employed_bees(ie);
    end

    % -4) calculate probability and apply roulette wheel selection scheme
    % to choose a food source for onlooker bee
    onlooker_bees();
    
    % I was rather seeing memorization here before scouting...
    
    % -8) if trail_i reaches max, abandon the source
    abandon=find(trail>=limit);
    for ia=abandon
        [x(:,ia),f(ia)]=scout_bees(ia);
        if debug
            disp([mfilename,':: site abandonment (',int2str(ia),')'])
            pause(0.25)
        end
    end
    
    memorize_best_food_source();
    
    if de_step
        differential_evolution_step()
    end
    if rem(obj.iter,verbose)==0 || obj.iter==1
        restart=1;
        disperse=dispersion(x,lb,ub);
        display_progress(restart,obj.iter,fbest,min(f),...
            disperse,obj.fcount,obj.optimizer);
    end
    stopflag=check_convergence(obj);
end

    function differential_evolution_step()
        [f,x,trail]=resort(f,x,trail);
        for ip=1:nparticles
            v=mutation(ip);
            u=crossover(v,ip);
            selection(u,ip)
        end
        function v=mutation(index)
            r=draw_indices(index,nparticles);
            x23=get_difference(x(:,r(2)),x(:,r(3)));
            v=x(:,r(1))+F*x23;
            v=reposition(v);
        end
        function u=crossover(v,index)
            u=x(:,index);
            crossers=rand(npar,1)<=cr;
            % ensure that at least one guy gets changed in case cr is too
            % small
            jrand=randperm(npar);
            jrand=jrand(1);
            crossers(jrand)=true;
            u(crossers)=v(crossers);
        end
    end
    function selection(u,index)
        fu=evaluate(u);
        improved=fu<f(index);
        if improved
            f(index)=fu;
            x(:,index)=u;
            trail(index)=0;
        else
            trail(index)=trail(index)+1;
        end
    end
    function ind=draw_indices(id,maxk)
        ind=[1:id-1,id+1:maxk];
        order=randperm(maxk-1);
        ind=ind(order);
    end
    function memorize_best_food_source()
        best=find(f==min(f),1,'first');
        if f(best)<fbest
            fbest=f(best);
            xbest=x(:,best);
        end
    end
    function v=reposition(v)
        bad=v<lb;
        v(bad)=lb(bad);
        bad=v>ub;
        v(bad)=ub(bad);
    end
    function employed_bees(index)
        u=rand(1,2);% <-- u=rand(npar,2);
        phi=-1+2*u;
        kk=draw_indices(index,n_empl);
        xi1=get_difference(x(:,index),x(:,kk(1)));
        v=x(:,index)+phi(:,1).*xi1;
        if abc>1
            xi2=get_difference(x(:,kk(2)),x(:,kk(3)));
            v=v+phi(:,2).*xi2;
        end
        v=reposition(v);
        selection(v,index);
    end
    function onlooker_bees()
        % -4) calculate probability and apply roulette wheel selection scheme
        % to choose a food source for onlooker bee
        fitness=f/(1+max(f)-min(f));
        fitness=1./(1+fitness);
        p=fitness/sum(fitness);
        cp=[0,cumsum(p)];
        cp(end)=1;
        for io=1:n_onlook
            isel=find(cp>=rand,1,'first')-1;
            employed_bees(isel);
        end      
    end
    function [xnew,fnew]=scout_bees(index)
        ublb=ub-lb;
        x1=x(:,index)+ublb/8.*randn(npar,1);
        x2=lb+ublb.*rand(npar,1);
        left=rand(npar,1);
        xnew=left.*x1+(1-left).*x2;
        fnew=evaluate(xnew);
    end
    function [fi]=evaluate(xi)
        fi=objective(xi);
        obj.fcount=obj.fcount+1;
    end
    function xab=get_difference(xa,xb)
%         d=distance(xa,xb);
%         if d<threshold
%             xb=lb+rand(npar,1).*(ub-lb);
%         end
        xab=xa-xb;
    end
end


function d=distance(xa,xb)
d=sqrt(sum(bsxfun(@minus,xa,xb).^2));
end