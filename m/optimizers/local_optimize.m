function [x,fx,obj]=local_optimize(Objective,x,fx,lb,ub,Options,varargin)
obj=struct('iterations',0,...
'MaxIter',inf,...
'start_time',clock,...
'MaxTime',inf,...
'funcCount',0,...
'MaxFunEvals',inf,...
'intial_step_size',1,...
'verbose',10,...
'optimizer',mfilename,...
'restart',0,...
'penalty',1e10,...
'threshold',1e-6);
Fields=fieldnames(obj);
for ii=1:size(Fields,1)
    if isfield(Options,Fields{ii,1})
        obj.(Fields{ii})=Options.(Fields{ii});
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
if isempty(fx)
    fx=evaluate(x);
end

u=zeros(size(x));
v = u;
par_id = 0;% -1;
vector_works = true;
vr = -obj.intial_step_size;
if isempty(fx)
    fx  = evaluate(x);
end
fxv = obj.penalty;
xv=x;
if obj.verbose %#ok<*BDSCI,*BDLGI>
    display_progress(obj.restart,obj.iterations,fx,fx,0,obj.funcCount,obj.optimizer)
end
if ~exist('ManualStoppingFile.txt','file')
    manual_stopping;
end
if abs(vr)<obj.threshold
    disp([mfilename,':: Peak found (iterations=',int2str(obj.iterations),',MaxIter=',int2str(obj.MaxIter),')'])
    obj.iterations=obj.MaxIter+1;
end
stopflag=check_convergence(obj);

while isempty(stopflag)
    obj.iterations=obj.iterations+1;
    if abs(vr)<2*obj.threshold
        inner_max_iter = 2*npar_v;
    else
        inner_max_iter =2;
    end
    
    inner_iter = 0;
    while fxv >= fx && inner_iter < inner_max_iter
        vector_works = false;
        if inner_iter == 0
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
        inner_iter=inner_iter+1;
    end
    
    if fxv >= fx
        fxv = obj.penalty;
        vr = .5*vr;
    else
        fx = fxv;
        x = xv;
        if inner_iter == 0
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
    if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
        display_progress(obj.restart,obj.iterations,fx,fx,vr,obj.funcCount,obj.optimizer)
    end
    if abs(vr)< obj.threshold
        disp([mfilename,':: Peak found (iterations=',int2str(obj.iterations),',MaxIter=',int2str(obj.MaxIter),')'])
        obj.iterations=obj.MaxIter+1;
    end
    stopflag=check_convergence(obj);
end

x=restore(x);

    function fy=evaluate(y)
        y=restore(y);
        fy=Objective(y,varargin{:});
        obj.funcCount=obj.funcCount+1;
    end

    function xx=restore(x)
        ncols=size(x,2);
        xx=nan(npar,ncols);
        xx(fixed,:)=fixed_vals(:,ones(ncols,1));
        xx(variable,:)=x;
    end
end