function [x,fx,fcount,outer_iter]=local_optimize_2(Objective,x,lb,ub,Options,varargin)
TolFun=[];intial_step_size=[];verbose=[];MaxIter=[]; MaxFunEvals=[];MaxTime=[];
optimizer=[];fcount=[];restart=[]; penalty=[];
Fields={'TolFun',1e-6
    'MaxIter',inf
    'MaxFunEvals',20000
    'MaxTime',3*60*60
    'intial_step_size',1
    'verbose',true
    'optimizer',mfilename
    'fcount',0
    'restart',0
    'penalty',1e10};
for ii=1:size(Fields,1)
    if isfield(Options,Fields{ii,1}) && ~isempty(Options.(Fields{ii,1}))
        eval([Fields{ii,1},'=Options.',Fields{ii,1},';'])
    else
        eval([Fields{ii,1},'=Fields{ii,2};'])
    end
end

npar=size(lb,1);
fixed=ub==lb;
fixed_vals=lb(fixed);
variable=find(~fixed);
npar_v=numel(variable);

% trim everything
lb=lb(variable);
ub=ub(variable);
if isempty(x)
    x=lb+(ub-lb).*rand(npar_v,1);
else
    x=x(variable);
end
fx=evaluate(x);

u=zeros(size(x));
v = u;
par_id = 0;% -1;
vector_works = true;
vr = -intial_step_size;
if isempty(fx)
    fx  = evaluate(x);
end
fxv = penalty;
xv=x;
outer_iter=0;
if verbose %#ok<*BDSCI,*BDLGI>
    display_progress(restart,outer_iter,fx,fx,0,fcount,optimizer)
end
t0=clock;
while abs(vr)>= TolFun && outer_iter<MaxIter && fcount < MaxFunEvals && etime(clock,t0)<MaxTime
    outer_iter=outer_iter+1;
    if abs(vr)<2*TolFun
        inner_max_iter = 2*npar_v;
    else
        inner_max_iter =2;
    end
    
    iter = 0;
    while fxv >= fx && iter < inner_max_iter
        vector_works = false;
        if iter == 0
            xv = x;
        else
            xv(par_id)= xv(par_id)-vr;
        end
        vr = -vr;
        
        if vr > 0
            par_id = par_id+1;
            if par_id>npar_v
                par_id=1;
            end % <---par_id = max(1,mod(par_id+1,npar_v));
        end
        xv(par_id) =xv(par_id)+ vr;
        xv(par_id)=recenter(xv(par_id),lb(par_id),ub(par_id));
        fxv  = evaluate(xv);
        iter=iter+1;
    end
    
    if fxv >= fx
        fxv = penalty;
        vr = .5*vr;
    else
        fx = fxv;
        x = xv;
        if iter == 0
            if vector_works
                u = u+v;
                v = v*2;
                xv = recenter(xv+v,lb,ub);
                vr = 2*vr;
            else
                u(par_id) = u(par_id)+vr;
                vr = 2*vr;
                xv(par_id) = xv(par_id)+vr;
                xv = recenter(xv,lb,ub);
            end
            fxv  = evaluate(xv);
        else
            xv = xv+u;
            xv(par_id) = xv(par_id) +vr;
            xv = recenter(xv,lb,ub);
            fxv  = evaluate(xv);
            if fxv >= fx
                u = 0*u; xv = x;
                u(par_id) = vr; vr= 2*vr;
                xv(par_id)= xv(par_id)+vr;
                xv = recenter(xv,lb,ub);
            else
                x= xv; fx = fxv;
                u(par_id) = u(par_id) +vr;
                v = 2*u; 
                vector_works = true;
                xv= recenter(xv+v,lb,ub);
                vr = sqrt(sum(v.^2));
            end
            fxv = evaluate(xv);
        end
    end
    if verbose
        display_progress(restart,outer_iter,fx,fx,vr,fcount,optimizer)
    end
end

x=restore(x);

    function fy=evaluate(y)
        y=restore(y);
        fy=Objective(y,varargin{:});
        fcount=fcount+1;
    end

    function xx=restore(x)
        ncols=size(x,2);
        xx=nan(npar,ncols);
        xx(fixed,:)=fixed_vals(:,ones(ncols,1));
        xx(variable,:)=x;
    end
end