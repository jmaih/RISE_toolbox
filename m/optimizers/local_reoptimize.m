function [xstore,fstore,fcount,outer_iter]=local_reoptimize(Objective,x,fx,lb,ub,Options,varargin)

threshold=[]; intial_step_size=[]; verbose=[]; max_iter=[];
optimizer=[]; fcount=[]; restart=[];

Fields={'threshold',1e-6
    'intial_step_size',1
    'verbose',50
    'max_iter',inf
    'optimizer',mfilename
    'fcount',0
    'restart',0};
for ii=1:size(Fields,1)
    if isfield(Options,Fields{ii,1})
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
if isempty(fx)
    fx=evaluate(x);
end

delta=intial_step_size;
max_delta=delta;
full_cycles=0;
vector=zeros(npar_v,1);
outer_iter=0;
fz=inf;
par_order=randperm(npar_v);
% par_order=1:npar_v;
par_count=0;
norm_delta=norm(delta);
fstore=[];
xstore=[];
while norm_delta>threshold && outer_iter<max_iter%#ok<*BDSCI>
    outer_iter=outer_iter+1;
    if norm_delta<2*threshold
        max_inner_iter = npar_v;
    else
        max_inner_iter =1;
    end
    inner_iter=0;
    while fz>=fx && inner_iter<max_inner_iter
        inner_iter=inner_iter+1;
        if inner_iter==1
            % reset z to a good start
            z=x;
        end
        par_count=par_count+1;
        if par_count>npar_v
            par_count=1;
            invalid=true;
            while invalid
                par_order=randperm(npar_v);
                if par_order(1)~=p_id
                    invalid=false;
                end
            end
            full_cycles=full_cycles+1;
            if isempty(fstore)||fx<min(fstore)
                fstore=[fstore,fx]; %#ok<AGROW>
                xstore=[xstore,x]; %#ok<AGROW>
            end
        end
        p_id=par_order(par_count);
        % try positive
        d=min(abs(delta),ub(p_id)-z(p_id));
        z(p_id)=z(p_id)+d;
        fz=evaluate(z);
        success=fz<fx;
        if ~success
            % try negative
            z(p_id)=z(p_id)-d;
            d=-min(abs(delta),z(p_id)-lb(p_id));
            z(p_id)=z(p_id)+d;
            fz=evaluate(z);
            success=fz<fx;
        end
        if success
            delta=d;
            if abs(delta)>abs(max_delta)
                max_delta=delta;
            end
        end
    end
    
    if fz>=fx % step too large?
        delta=.5*delta;
        fz=inf;
    else
        % step profitable. exploit more
        fx=fz;
        x=z;
        
        % try the full vector
        vector_success=vector_orientation(fx);
        
        if ~vector_success
            % further try the successful parameter
            individual_orientation(fx);
        end
    end
    norm_delta=norm(delta);
    if ismember(rem(outer_iter,verbose),[0,1])
        display_progress(restart,outer_iter,fx,fx,norm_delta,fcount,optimizer)
    end
end
display(max_delta)
display(full_cycles);

[fstore,xstore]=resort(fstore,restore(xstore));


    function individual_orientation(f_old)
        % a small step didn't work, try a bigger one and change
        % direction
        vector=0*vector;
        z=x;
        vector(p_id)=delta;
        delta=2*delta;
        z(p_id)=z(p_id)+delta;
        fz=evaluate(z);
        promising=fz<f_old;
        zstar=z;
        fstar=fz;
        while promising
            % the big step worked. try an even bigger one
            vector(p_id)=vector(p_id)+delta;
            delta=2*delta;
            z(p_id)=z(p_id)+delta;
            z(p_id)=recenter(z(p_id),lb(p_id),ub(p_id));
            fz=evaluate(z);
            promising=fz<f_old;
            if promising
                zstar=z;
                fstar=fz;
            end
            f_old=fz;
        end
        if abs(delta)>abs(max_delta)
            max_delta=delta;
        end
        z=zstar;
        fz=fstar;
    end

    function success=vector_orientation(f_old)
        % try whole vector
        z=z+vector;
        z(p_id)=z(p_id)+delta;
        z=recenter(z,lb,ub);
        fz=evaluate(z);
        success=false;
        if fz<fx
            success=true;
            % the vector is promising
            vector(p_id)=vector(p_id)+delta;
            v=2*vector;
            z=z+v;
            z=recenter(z,lb,ub);
            fz=evaluate(z);
            promising=fz<f_old;
            zstar=z;
            fstar=fz;
            vstar=v;
            while promising
                vector=vector+v;
                v=2*v;
                z=z+v;
                z=recenter(z,lb,ub);
                fz=evaluate(z);
                promising=fz<f_old;
                if promising
                    zstar=z;
                    fstar=fz;
                    vstar=v;
                end
                f_old=fz;
            end
            delta=norm(vstar);
            if delta>abs(max_delta)
                max_delta=delta;
            end
            z=zstar;
            fz=fstar;
        end
    end

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
