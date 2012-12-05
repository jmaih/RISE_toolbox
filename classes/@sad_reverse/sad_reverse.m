classdef sad_reverse
    properties
        name
        args
    end
    properties(SetAccess = protected)
        ready ={'plus','times','divide','minus','power'}
        % holds the list of signs such that if used on a specific object
        % parentheses are not required
    end
    methods
        function obj=sad_reverse(name,arg)
            if nargin
                if isa(name,'sad_reverse')
                    obj=name;
                else
                    if nargin<2, arg=[]; end
                    if isempty(arg),arg={}; end
                    if ~iscell(arg),arg={arg}; end
                    obj.name=name;
                    % algebraic simplification at node construction
                    obj.args=arg;
                    if ~isa(obj.name,'double')
                        % use the functional forms of ge gt le lt ne
                        if ismember(obj.name,{'ldivide','mldivide',...
                                'mrdivide','rdivide','mtimes','times',...
                                'power','mpower'})
                            obj.ready={'plus','times','divide','minus'};
                        elseif ismember(obj.name,{'minus','plus'})
                            obj.ready={'plus'};
                        end
                    end
                end
            end
        end
        function obj=plus(u,v)    ,obj=sad_reverse('plus',{u,v}); end
        function obj=uplus(u)     ,obj=sad_reverse('uplus',{u}); end
        function obj=minus(u,v)   ,obj=sad_reverse('minus',{u,v}); end
        function obj=uminus(u)    ,obj=sad_reverse('uminus',{u}); end
        function obj=times(u,v)   ,obj=sad_reverse('times',{u,v}); end
        function obj=mtimes(u,v)  ,obj=sad_reverse('mtimes',{u,v}); end
        function obj=power(u,v)   ,obj=sad_reverse('power',{u,v}); end
        function obj=mpower(u,v)  ,obj=sad_reverse('mpower',{u,v}); end
        function obj=rdivide(u,v) ,obj=sad_reverse('rdivide',{u,v}); end
        function obj=mrdivide(u,v),obj=sad_reverse('mrdivide',{u,v}); end
        function obj=ldivide(u,v) ,obj=sad_reverse('ldivide',{u,v}); end
        function obj=mldivide(u,v),obj=sad_reverse('mldivide',{u,v}); end
        function obj=exp(u)       ,obj=sad_reverse('exp',{u}); end
        function obj=log(u)       ,obj=sad_reverse('log',{u}); end
        function obj=log10(u)     ,obj=sad_reverse('log10',{u}); end
        function obj=cos(u)       ,obj=sad_reverse('cos',{u}); end
        function obj=acos(u)      ,obj=sad_reverse('acos',{u}); end
        function obj=cosh(u)      ,obj=sad_reverse('cosh',{u}); end
        function obj=sin(u)       ,obj=sad_reverse('sin',{u}); end
        function obj=asin(u)      ,obj=sad_reverse('asin',{u}); end
        function obj=sinh(u)      ,obj=sad_reverse('sinh',{u}); end
        function obj=tan(u)       ,obj=sad_reverse('tan',{u}); end
        function obj=atan(u)      ,obj=sad_reverse('atan',{u}); end
        function obj=tanh(u)      ,obj=sad_reverse('tanh',{u}); end
        function obj=min(u,v)     ,obj=sad_reverse('min',{u,v}); end
        function obj=max(u,v)     ,obj=sad_reverse('max',{u,v}); end
        function obj=sum(u,v)     ,obj=sad_reverse('sum',{u,v}); end
        function obj=normpdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=sad_reverse('normpdf',{u,v},w);
        end
        function obj=normcdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=sad_reverse('normcdf',{u,v},w);
        end
        function obj=abs(u)       ,obj=sad_reverse('abs',{u}); end
        function obj=isreal(u)    ,obj=sad_reverse('isreal',{u}); end
        function obj=sqrt(u)      ,obj=sad_reverse('sqrt',{u}); end
        function obj=norm(u)      ,obj=sad_reverse('norm',{u}); end
        function obj=gt(u,v)      ,obj=sad_reverse('gt',{u,v}); end
        function obj=ge(u,v)      ,obj=sad_reverse('ge',{u,v}); end
        function obj=lt(u,v)      ,obj=sad_reverse('lt',{u,v}); end
        function obj=le(u,v)      ,obj=sad_reverse('le',{u,v}); end
        function obj=sign(u)      ,obj=sad_reverse('sign',{u}); end
        function string=char(this)
            n=numel(this);
            for i0=1:n
                obj=this(i0);
                % char itself is already taken care of
                if isa(obj.name,'double')
                    thisString=sprintf('%0.10g',obj.name); % <-- thisString=num2str(obj.name,10);
                elseif isempty(obj.args) % variable
                    thisString=obj.name; % this should be a char
                else
                    func_name=obj.name;
                    if ismember(func_name,{'mldivide','mrdivide','mpower','mtimes'})
                        func_name=func_name(2:end);
                    end
                    args_=reprocess_arguments(obj.args);
                    switch func_name
                        case {'normpdf','normcdf'}
                            thisString=strcat(func_name,'(',args_{1},',',args_{2},',',args_{3},')');
                        case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                                'tan','atan','tanh','abs','sqrt','isreal','sign'}
                            thisString=strcat(func_name,'(',args_{1},')');
                        case {'min','max','gt','lt','ge','le'}
                            thisString=strcat(func_name,'(',args_{1},',',args_{2},')');
                        case {'uminus','uplus'}
                            thisString=rise_algebra_cat(func_name,args_{1});
                        case {'plus','minus','times','power','rdivide'}
                            thisString=rise_algebra_cat(func_name,args_{1},args_{2});
                        otherwise
                            error([func_name,' is undefined for objects of class ',mfilename])
                    end
                end
                if n>1
                    if i0==1
                        cellmode=iscell(thisString);
                        if cellmode
                            string=cell(n,numel(thisString));
                        else
                            string=cell(size(this));
                        end
                    end
                    if cellmode
                        string(i0,:)=thisString;
                    else
                        string{i0}=thisString;
                    end
                else
                    string=thisString;
                end
            end
            function args=reprocess_arguments(args)
                for ii=1:numel(args)
                    if isnumeric(args{ii})
                        args{ii}=sprintf('%0.10g',args{ii}); % <-- args{ii}=num2str(args{ii},10);
                    else
                        check=is_ready(args{ii});
                        args{ii}=char(args{ii});
                        if ~check
                            args{ii}=strcat('(',args{ii},')');
                        end
                    end
                end
                function flag=is_ready(arg)
                    alien=~ismember(func_name,{'minus','times','power','rdivide','ldivide'});
                    flag=alien||...
                        (strcmp(func_name,'minus') && ii==1)||...
                        ismember(func_name,arg(1).ready);
                    % there may be many arguments but they,ll all share the
%                     % same functor 
                end
            end
        end
        function print(tree)
            for it=1:numel(tree)
                if ~isempty(tree(it).ref)
                    fprintf(1,'%s\n',[tree(it).ref,' ---> ',tree(it).name]);
                end
                tree(it).args=reprocess_arguments(tree(it).args);
                for iarg=1:numel(tree(it).args)
                    if ~ischar(tree(it).args{iarg})
                        print(tree(it).args{iarg})
                    end
                end
            end
        end
        function this=diff(obj,wrt)
            % check that I do not need to output the object thanks to the
            % handle property
            if iscell(wrt)
                wrt=[wrt{:}];
            end
            if ~isa(wrt,'sad_reverse')
                error([mfilename,':: second argument must be sad_reverse object'])
            end
            wrt=wrt(:)';
            this=mydiff(obj,wrt);
        end
        function obj=repmat(this,M,N)
            if nargin == 2
                if isscalar(M)
                    siz = [M M];
                else
                    siz = M;
                end
            else
                siz = [M N];
            end
            indexes=cell(1,numel(siz));
            for ind=1:numel(siz)
                indexes{ind}=1:siz(ind);
            end
            obj=sad_reverse.empty(0);
            obj(indexes{:})=this;
        end
    end
    methods(Access=private)
        function this=mydiff(obj,wrt)
            nobj=numel(obj);
            this=sad_reverse.empty(0);
            guizo=sad_reverse(0);
            for iobj=1:nobj
                myobj=obj(iobj);
                if isnumeric(myobj.name)
                    tmp=guizo;
                    % derivatives =0 as above
                elseif isempty(myobj.args)
                    nwrt=numel(wrt);
                    tmp=repmat(guizo,1,nwrt);
                    wrt_names={wrt.name};
                    tmp(strcmp(myobj.name,wrt_names))=sad_reverse(1);
                else
                    tmp=differentiation_engine(myobj,wrt);
                end
                this(iobj,:)=tmp;
            end
        end
    end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
    end
end

function this=differentiation_engine(myobj,wrt)
args_=myobj.args;
for iarg=1:numel(args_)
    args_{iarg}=sad_reverse(args_{iarg});
end
u=args_{1};
du=mydiff(u,wrt);
switch myobj.name
    case {'gt','ge','lt','le','sign'}
        this=sad_reverse(0);
    case 'plus'
        v=args_{2};
        dv=mydiff(v,wrt);
        this=du+dv;
    case 'uplus'
        this=du;
    case 'minus'
        v=args_{2};
        dv=mydiff(v,wrt);
        this=du-dv;
    case 'uminus'
        this=-du;
    case {'mtimes','times'}
        v=args_{2};
        dv=mydiff(v,wrt);
        upv=du*v;
        vpu=dv*u;
        this=upv+vpu;
    case {'power','mpower'}
        v=args_{2};
        this=v*du*u^(v-1);
        if ~isnumeric(v.name)
            dv=mydiff(v,wrt);
            this=this+dv*log(u)*u^v;
        end
    case {'rdivide','mrdivide'}
        v=args_{2};
        dv=mydiff(v,wrt);
        this=(du*v-dv*u)/v^2; % only when v~=0
    case {'ldivide','mldivide'}
        v=args_{2};
        dv=mydiff(v,wrt);
        this=(u*dv-v*du)/u^2; % only when u~=0
    case 'exp'
        this=du*myobj;
    case 'log'
        this=du/u;
    case 'log10'
        this=(du/u)/log(10);
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
    case 'min'
        v=args_{2};
        muv=u<v;
        this=muv*u+(1-muv)*v;
    case 'max'
        v=args_{2};
        muv=u>v;
        this=muv*u+(1-muv)*v;
    case 'sum'
        this=sum(du);
    case 'normpdf'
        v=args_{2};
        w=args_{3};
        this=-du/w*(u-v)/w*myobj;
    case 'normcdf'
        v=args_{2};
        w=args_{3};
        this=du*normpdf(u,v,w);
    case 'abs'
        this=du*(-u<0+u>0);
    case 'isreal'
        this=isreal(u)*du;
    case 'sqrt'
        this=du/(2*sqrt(myobj));
    case 'norm' % this would not work!
        this=sum(u.*du)/norm(u);
    otherwise
        error([myobj.name,' is unknown type of operator'])
end

end

% function algebraic_simplification_at_node_construction()
% if is_node_construction
%     % use the expression optimizer to simply the expression...
% numerics=cellfun(@isnumeric,{args.name});
% if any()
% a ? 1 ? a 
% a ? ?1 ? ?a
% a ? 0 ? 0 
% a±0 ? a
% a/a ? 1 
% a/?1 ? ?a
% a?a ? 0 
% f(c0) ? Constant(f(c0))
% c0 ?-+/ c1 ? Constant(c0 ? c1)
% c0 ±c1 ? Constant(c0 ±c1)
% end
