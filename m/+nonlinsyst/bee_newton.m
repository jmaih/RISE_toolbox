function [x,fval,exitflag] = bee_newton(f,x0,opts,lb,ub,varargin)

if nargin<3
    
    opts=[];
    
end

fcount=0;

nv=length(x0);

if isempty(lb)
    
    lb=-10000*ones(nv,1);
    
    ub=-lb;
    
end

nw=20;

opts=nonlinsyst.set_defaults(opts);

is_verbose=strcmp(opts.Display,'iter');

exitflag=1;

xf0=struct('x',nan,'fval',inf,'fnorm',inf);

x0=encode(x0);

if is_solution(x0)
    
    x=x0.x;
    
    fval=x0.fval;
    
    return
    
end

% try a quick Newton just in case the problem is simple
%------------------------------------------------------
[x,fval,exitflag] = nonlinsyst.mynewton(@objective,x0.x,opts);

if exitflag==1
    
    return
    
end

% set random number generator
scurr=rng();

rng(1942)

iter=0;

MaxIter=50;

newtOpts=opts;

newtOpts.Display='none';

% start at the best value for the previous newton
%-------------------------------------------------
x0.x=x; x0.fval=fval; x0=encode(x0);

while iter<MaxIter
    
    iter=iter+1;
    
    x0=bee_phase(x0,nw,lb,ub,@encode,xf0,newtOpts);
    
    if is_verbose
        
        fnormBees=x0(1).fnorm;
        
    end
    
    [x0,exitflag]=newton_phase(x0);
    
    if is_verbose
        
        fnormNewton=x0(1).fnorm;
        
        fprintf(1,'iter %d  fnorm(bees)= %8.4g  fnorm(newton)= %8.4g\n',...
            iter,fnormBees,fnormNewton);
        
    end
    
    [x0,order]=re_order(x0);
    
    exitflag=exitflag(order);
    
    if any(exitflag(1)==1)
        
        break
        
    end
           
end

% set solution
x=x0(1).x;

fval=x0(1).fval;

% reset random number generator
rng(scurr)


    function [x0,exitflag]=newton_phase(x0)
        
        exitflag=nan(1,nw);
        
        for iv=1:nw
            
            [x0(iv).x,x0(iv).fval,exitflag(iv)] = nonlinsyst.mynewton(@objective,x0(iv).x,newtOpts);
            
            x0(iv)=encode(x0(iv));
            
            if exitflag(iv)==1
                
                break
                
            end
            
        end
        
    end

    function e=encode(x)
        
        if isempty(x)
            
            e=xf0(1:0);
            
        elseif isnumeric(x)
            
            e=xf0;
            
            e.x=x;
            
            e.fval=objective(x);
            
            e.fnorm=nonlinsyst.norm(e.fval);
            
        elseif isstruct(x)
            
            e=x;
            % recalculate norm
            e.fnorm=nonlinsyst.norm(e.fval);
            
        end
        
        if ~is_feasible(e)
            
            e.fnorm=inf;
            
        end
        
    end

    function flag=is_solution(xf)
        
        flag=xf.fnorm<=opts.TolFun;
        
    end

    function o=objective(x)
        
        o=f(x,varargin{:});
        
        fcount=fcount+1;
        
    end

end

function flag=is_feasible(xf)

flag=isfinite(xf.fnorm);

end


function [xf,order]=re_order(xf)

[~,order]=sort([xf.fnorm]);

xf=xf(order);

end


function  xf=bee_phase(x0,nw,lb,ub,encode,xf0,opts)

discard_duplicates()

is_verbose=strcmp(opts.Display,'iter');

MaxIter=100;

max_rounds=25;

nv=numel(x0(1).x);

food_number=round(.5*nw);

max_genes_change = 1; % number of genes change

[xf,trials]=set_warriors(x0,nw);

xbest=xf(1);

iter=0;

max_trials=min(round(0.5*MaxIter),50);

% Print header.
if is_verbose
    
    fprintf('\n Iteration     norm f(x)\n');
    
    formatstr = '%5.0f   %15.6g\n';
    
end

stopflag=check_convergence();

while ~stopflag
    
    iter=iter+1;
        
    send_employed_bees()
    
    send_onlooker_bees()
    
    memorize_best_source()
    
    send_scout_bees()
    
    display_progress()
    
    stopflag=check_convergence();
    
end

fnorm=[xf.fnorm];

if ~any(fnorm-xbest.fnorm==0)
    
    xf=re_order([xf,xbest]);
    
    xf=xf(1:nw);
    
end

    function discard_duplicates()
        
        sqrtEps=opts.TolX;
        
        n0=numel(x0);
        
        if n0==1
            
            return
            
        end
        
        xx=[x0.x];
        
        xfnorm=num2cell([x0.fnorm]);
        
        discard=false(1,n0);
        
        for icol=2:n0
            
            if discard(icol-1)
                
                continue
                
            end
            
            ddd=max(abs(bsxfun(@minus,xx(:,icol-1),xx(:,icol:end))));
            
            bad=find(ddd<sqrtEps);
            
            discard(bad+icol-1)=true;
            
        end
        
        x0=x0(~discard);
        
    end


    function bees_reorder()
        
        [xf,order]=re_order(xf);
        
        trials=trials(order);
        
    end


    function display_progress()
        
        if ~is_verbose
            
            return
            
        end
        
        fprintf(1,formatstr,iter,xbest.fnorm);
        
    end


    function flag=check_convergence()
                
        flag=iter>=MaxIter;
            
    end


    function send_scout_bees()
        
        renew=find(trials>=max_trials);
        
        if isempty(renew)
            
            return
            
        end
        
        [xf(renew),trials(renew)]=set_warriors(x0(1:0),numel(renew));
        
    end


    function memorize_best_source()
        
        bees_reorder()
        
        if xf(1).fnorm < xbest.fnorm
            
            xbest=xf(1);
            
        end
        
    end


    function send_onlooker_bees()
        
        fitness=utils.optim.compute_fitness([xf.fnorm]);
        
        prob=(0.9.*fitness./max(fitness))+0.1;
        
        t=0;
        
        ii=1;
        
        while t<food_number
            
            if rand<prob(ii)
                
                t=t+1;
                
                generate_mutant(ii);
                
            end
            
            ii=ii+1;
            
            if ii==food_number+1
                
                ii=1;
                
            end
            
        end
        
    end


    function send_employed_bees()
        
        for ii=1:food_number
            
            generate_mutant(ii)
            
        end
        
        bees_reorder()
        
    end


    function generate_mutant(ii)
        % generate a solution for index ii
        % a randomly chosen solution different from ii is used for producing a mutant
        while 1
            
            donor_id=min(food_number,fix(rand*food_number)+1);
            
            if donor_id~=ii
                
                break
                
            end
            
        end
        
        xmut=xf(ii).x;
        
        %  pick the parameters to change in the new solution
        nch=min(randi(ceil(0.5*nv)),max_genes_change);
        
        change=min(fix(rand(nch,1)*nv)+1,nv);
        
        xmut(change)=xmut(change)+(xmut(change)-xf(donor_id).x(change))*2.*(rand(nch,1)-.5);
        
        xmut(change)=utils.optim.recenter(xmut(change),lb(change),ub(change));
        
        mutant=encode(xmut);
        
        if mutant.fnorm<xf(ii).fnorm
            
            xf(ii)=mutant;
            
            trials(ii)=0;
            
        else
            
            trials(ii)=trials(ii)+1;
            
        end
        
    end


    function [xf,trials]=set_warriors(x0,nw)
        
        trials=zeros(1,nw);
        
        n0=numel(x0);
        
        feas=false(1,n0);
        
        for iw=1:n0
            
            feas(iw)=is_feasible(x0(iw));
            
        end
        
        patch_lb=max(lb,-1);
        patch_ub=min(ub,1);
        
        xf=x0(feas);
        
        istart=numel(xf)+1;
        
        for iw=istart:nw
            
            xfd=xf0;
            
            iround=0;
            
            while ~is_feasible(xfd) && iround<max_rounds
                
                xfd=encode(patch_lb+(patch_ub-patch_lb).*rand(nv,1));
                % xfd=encode(lb+(ub-lb).*rand(nv,1));
                
                iround=iround+1;
                
            end
            
            if iround==max_rounds
                
                error('failed to find warriors')
                
            end
            
            xf(iw)=xfd;
            
        end
        
    end


end