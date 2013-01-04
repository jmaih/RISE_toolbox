function this=diff(obj,wrt)
% check that I do not need to output the object thanks to the
% handle property
if iscell(wrt)
    wrt=[wrt{:}];
end
if ~isa(wrt,'sadiff')
    error([mfilename,':: second argument must be sadiff object'])
end
wrt=wrt(:)';
zeroObj=sadiff(0);
this=mydiff(obj,wrt,zeroObj);
% this=expand(this);
end

function this=mydiff(obj,wrt,zeroObj)
if nargin<3
    zeroObj=sadiff(0);
end
nobj=numel(obj);
this=sadiff.empty(0);
nwrt=numel(wrt);
wrt_names={wrt.func};
for iobj=1:nobj
    myobj=obj(iobj);
    if isnumeric(myobj.func)
        tmp=zeroObj;%(ones([1,nwrt]));
        % derivatives =0 as above
    elseif isempty(myobj.args)
        loc= strcmp(myobj.func,wrt_names);
        tmp=zeroObj(ones([1,nwrt])); %<--repmat(zeroObj,1,nwrt);
        % then this variable is differentiable
        if any(loc)
        thisOrder=wrt(loc).order;
        tmp(thisOrder)=sadiff(1);
        end
    else
        tmp=differentiation_engine(myobj,wrt,zeroObj);
    end
    this(iobj,:)=tmp;
end
end

function this=differentiation_engine(myobj,wrt,zeroObj)
args_=myobj.args;
nargs=numel(args_);
for iarg=1:numel(args_)
    args_{iarg}=sadiff(args_{iarg});
end
u=args_{1};
du=mydiff(u,wrt,zeroObj);
if nargs>1
    v=args_{2};
        dv=mydiff(v,wrt,zeroObj);
    if nargs>2
        w=args_{3};
        %                     dw=mydiff(w,wrt);
        if nargs>3
            error('functions with more that 3 arguments not implemented')
        end
    end
end

switch myobj.func
    case '+'
        this=du+dv;
    case 'uplus'
        this=du;
    case '-'
        this=du-dv;
    case 'uminus'
        this=-du;
    case {'*'}
        upv=du*v;
        vpu=dv*u;
        this=upv+vpu;
    case {'^'}
        this=v*du*u^(v-1);
        if ~isnumeric(v.func)
            this=this+dv*log(u)*u^v;
        end
    case {'/'}
        this=(du*v-dv*u)/v^2; % only when v~=0
    case 'exp'
        this=du*myobj;
    case 'log'
        this=du/u;
    case 'cos'
        this=-du*sin(u);
    case 'acos'
        this=-du/sqrt(1-u^2);
    case 'cosh'
        this=du*sinh(u);
    case 'sin'
        this=du*cos(u);
    case 'asin'
        this=du/sqrt(1-u^2);
    case 'sinh'
        this=du*cosh(u);
    case 'tan'
        this=du/(cos(u))^2;
    case 'atan'
        this=du/(1+u^2);
    case 'tanh'
        this=du/(1-u^2);
    case 'normpdf'
        this=-du/w*(u-v)/w*myobj;
    case 'normcdf'
        this=du*normpdf(u,v,w);
    case 'sqrt'
        this=du/(2*sqrt(myobj));
    otherwise
        error([myobj.func,' is unknown type of operator'])
end

end
