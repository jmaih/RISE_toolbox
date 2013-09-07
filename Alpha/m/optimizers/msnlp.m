function [xbest,fbest,exitflag,output]=msnlp(objective,x0,lb,ub,options,varargin)

default_options = struct('dwait',20,'mwait',20,'dfuser',0.8,'mfuser',0.2,...
    'n1',200,'n2',1000,'MaxIter',150000,'MaxFunEvals',2000000,...
    'local_optimizer',@fmincon,'TolX',1e-4,'TolFun',1e-4,...
    'Display','iter','nonlcon',[]);

% If just 'defaults' passed in, return the default options in X
if nargin==0
    xbest=default_options;
    return
end
local_options={'dwait','mwait','dfuser','mfuser','n1','n2','local_optimizer','nonlcon'};

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

dwait=default_options.dwait;
mwait=default_options.mwait;
dfuser=default_options.dfuser;
mfuser=default_options.mfuser;
n1=default_options.n1;
n2=default_options.n2;
local_optimizer=default_options.local_optimizer;
if isa(local_optimizer,'function_handle')
    optimizer_name=func2str(local_optimizer);
else
    optimizer_name=local_optimizer;
    local_optimizer=str2func(local_optimizer);
end
nonlcon=default_options.nonlcon;

basin_radius=[];
drejctr=[];
dftemp=[];
mftemp=inf;
nlocals=0;
funcCount=0;
mrejctr=0;

npar=numel(x0);
x0=wrapper(x0);

xf_local=wrapper();
% stage 1

optim_opt=rmfield(default_options,local_options);
optim_opt.MaxFunEvals=ceil(optim_opt.MaxFunEvals/(n2+1));
optim_opt.MaxIter=ceil(optim_opt.MaxIter/(n2+1));
matlab_flag=ismember(optimizer_name,{'fmincon','fminsearch','fminunc'});
if matlab_flag
    Problem=struct('objective',@optimization_wrapper,...
        'x0',x0.x,...
        'lb',lb,...
        'ub',ub,...
        'A',[],...
        'B',[],...
        'Aeq',[],...
        'Beq',[],...
        'solver',optimizer_name,...
        'nonlcon',nonlcon,...
        'options',optim_opt);
    xf=local_optimizer(Problem);
else
    xf=local_optimizer(@optimization_wrapper,x0.x,optim_opt);
end
xf=wrapper(xf);
update_locals(x0,xf)

xstar=[];
for i=1:n1
    xi=next_candidate();
    if isempty(xstar)
        xstar=xi;
    else
        xstar=selection_process(xstar,xi);
    end
end

if matlab_flag
    Problem.x0=xstar.x;
    xf=local_optimizer(Problem);
else
    xf=local_optimizer(@optimization_wrapper,xstar.x,optim_opt);
end
xf=wrapper(xf);

update_locals(xstar,xf)

threshold=xstar.f;

% stage 2

for i=1:n2
    xi=next_candidate();
    dstatus=distance_filter(xi) ;
    if dstatus
        mstatus=merit_filter(xi);
        if mstatus
            if matlab_flag
                Problem.x0=xi.x;
                xf=local_optimizer(Problem);
            else
                xf=local_optimizer(@optimization_wrapper,xi.x,optim_opt);
            end
            xf=wrapper(xf);
            update_locals(xi,xf);
        end
    end
end

xf_local=sort_population(xf_local);
xbest=xf_local(1).x;
fbest=xf_local(1).f;
exitflag=1;
output=struct('funcCount',funcCount,'iterations',n1+n2);

%--------------------------------------------------------------------------
    function update_locals(xf_start,xf_final)
        % process and stor solver output and produce penalty weights w
        if isfeasible(xf_final)
            d=distance(xf_start.x,xf_final.x);
            [flag,loc]=isnewlocal(xf_final);
            if flag
                xf_local=[xf_local,xf_final];
                basin_radius=[basin_radius,d];
                nlocals=nlocals+1;
                dftemp=[dftemp,0];
                drejctr=[drejctr,0];
            else
                basin_radius(loc)=max(d,basin_radius(loc));
            end
            % basin overlap exclusion
            for ii=1:nlocals
                for jj=1:nlocals
                    if jj==ii
                        continue
                    end
                    dij=distance(xf_local(ii).x,xf_local(jj).x);
                    if dij<basin_radius(ii)+basin_radius(jj)
                        f=dij/(basin_radius(ii)+basin_radius(jj));
                        basin_radius(ii)=f*basin_radius(ii);
                        basin_radius(jj)=f*basin_radius(jj);
                        drejctr([ii,jj])=0;
                    end
                end
            end
        end
    end
%--------------------------------------------------------------------------
    function dstatus=distance_filter(xnew)
        dstatus=true;
        for iloc=1:nlocals
            if drejctr(iloc)==dwait
                df=min(dfuser,dftemp(iloc));
                basin_radius(iloc)=df*basin_radius(iloc);
                drejctr(iloc)=0;
                dftemp(iloc)=0;
            end
            di=distance(xnew.x,xf_local(iloc).x);
            if di<dfuser*basin_radius(iloc)
                dstatus=false;
                drejctr(iloc)=drejctr(iloc)+1;
                dftemp(iloc)=max([dftemp(iloc),di/basin_radius(iloc),0.5]);
            end
        end
    end
%--------------------------------------------------------------------------
    function mstatus=merit_filter(xnew)
        mstatus=false;
        if xnew.f<threshold
            mstatus=true;
            threshold=xnew.f;
            mrejctr=0;
            mftemp=inf;
        else
            mrejctr=mrejctr+1;
            mftemp=min(mftemp,(xnew.f-threshold)/(1+abs(threshold)));
            if mrejctr==mwait
                mrejctr=0;
                mf=max(mfuser,mftemp);
                threshold=threshold+mf*(1+abs(threshold));
            end
        end
    end

%--------------------------------------------------------------------------
    function flag=isfeasible(xf_final)
        flag=xf_final.violstrength==0;
    end
%--------------------------------------------------------------------------
    function [flag,loc]=isnewlocal(individual)
        loc=find(abs([xf_local.f]-individual.f)<1e-6,1,'first');
        flag=isempty(loc);
    end
%--------------------------------------------------------------------------
    function f=optimization_wrapper(x)
        this=evaluate_individual(x,objective,lb,ub,nonlcon,varargin{:});
        f=this.f;
    end
%--------------------------------------------------------------------------
    function this=wrapper(bird)
        if nargin==0
            this=evaluate_individual();
        else
            this=evaluate_individual(bird,objective,lb,ub,nonlcon,varargin{:});
            funcCount=funcCount+1;
        end
    end
    function c=next_candidate()
        c=wrapper(lb+(ub-lb).*rand(npar,1));
    end
end