classdef planar
    properties
        func
        args
        incidence % tells which variable appears
        is_vector=false
        lineage={} % will serve to trace back all generations of derivatives
        prototype
        location % will hold position (column) where the derivative will be stored
    end
    methods
        % constructor
        %------------
        function obj=planar(f,a)
            if nargin>0
                if strcmp(class(f),'planar') %#ok<*STISA>
                    obj=f;
                else
                    obj.func=f;
                    if nargin>1
                        if ~iscell(a)
                            a={a};
                        end
                        obj.args={a};
                    end
                    obj.is_vector=isnumeric(f) && numel(f)>1;
                end
            end
        end
        % overloaded operators
        %---------------------
        function obj=abs(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=abs(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='abs';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=acos(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=acos(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='acos';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=acosh(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=acosh(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='acosh';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=and(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=and(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='and';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=asin(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=asin(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='asin';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=asinh(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=asinh(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='asinh';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=atan(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=atan(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='atan';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=atanh(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=atanh(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='atanh';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=cos(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=cos(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='cos';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=cosh(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=cosh(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='cosh';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=cot(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=cot(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='cot';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=eq(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=eq(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='eq';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=erf(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=erf(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='erf';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=exp(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=exp(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='exp';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=ge(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=ge(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='ge';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=gt(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=gt(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='gt';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=if_elseif(varargin)
            obj=[];
            n=length(varargin);
            constants=cell(1,n);
            nconst=0;
            guy_is_planar=false(1,n);
            for iarg=1:n
                if strcmp(class(varargin{iarg}),'planar')
                    guy_is_planar(iarg)=true;
                    if isempty(obj)
                        obj=varargin{iarg}.prototype;obj.prototype=obj;
                    end
                    if isnumeric(varargin{iarg}.func)
                        nconst=nconst+1;
                        constants{iarg}=varargin{iarg}.func;
                    end
                else
                    nconst=nconst+1;
                    constants{iarg}=varargin{iarg};
                end
            end
            for iarg=find(~guy_is_planar)
                tmp=varargin{iarg};
                varargin{iarg}=obj; varargin{iarg}.func=tmp;
            end
            if nconst==n
                obj.func=if_elseif(constants{:});
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='if_elseif';
                obj.args=varargin;
                for iarg=1:n
                    obj.is_vector=obj.is_vector||varargin{iarg}.is_vector;
                    if ~isempty(varargin{iarg}.incidence)
                        if isempty(obj.incidence)
                            obj.incidence=varargin{iarg}.incidence;
                        else
                            obj.incidence=obj.incidence|varargin{iarg}.incidence;
                        end
                    end
                end
            end
        end
        function obj=if_then_else(a,b,c)
            obj=[];
            if strcmp(class(a),'planar'),obj=a.prototype;obj.prototype=obj; end
            if strcmp(class(b),'planar') && isempty(obj),obj=b.prototype;obj.prototype=obj; end
            if strcmp(class(c),'planar') && isempty(obj),obj=c.prototype;obj.prototype=obj; end
            if ~strcmp(class(a),'planar'),atmp=a; a=obj;a.func=atmp; end
            if ~strcmp(class(b),'planar'),btmp=b; b=obj;b.func=btmp; end
            if ~strcmp(class(c),'planar'),ctmp=c; c=obj;c.func=ctmp; end
            if isnumeric(a.func) && isnumeric(b.func) && isnumeric(c.func)
                obj.func=if_then_else(a.func,b.func,c.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='if_then_else';
                obj.args={a,b,c};
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
                if ~isempty(c.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=c.incidence;
                    else
                        obj.incidence=obj.incidence|c.incidence;
                    end
                end
                obj.is_vector=a.is_vector||b.is_vector||c.is_vector;
            end
        end
        function obj=le(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=le(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='le';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=log(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=log(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='log';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=log10(a)
            obj=log(a)/log(10);
        end
        function obj=lt(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=lt(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='lt';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=max(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=max(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='max';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=min(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=min(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='min';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=minus(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=minus(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            elseif is_zero(b)
                % x - 0 = x
                obj=a;
            elseif is_zero(a)
                % 0 - x = - x
                obj=-b;
            elseif ischar(b.func) && strcmp(b.func,'uminus')
                % x - (-y) = x + y
                obj=a+b.args{1};
            else
                obj.func='minus';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=mpower(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=power(a.func,b.func); % note the vectorized form
                obj.is_vector=numel(obj.func)>1;
            elseif is_zero(b)
                % x^0 = 1
                obj.func=1;
                %	obj.is_vector=numel(obj.func)>1;
            elseif is_one(b)
                % x^1 = x
                obj=a;
            elseif is_one(a)
                % 1^x = 1
                obj.func=1;
                %	obj.is_vector=numel(obj.func)>1;
            else
                obj.func='mpower';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=mrdivide(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=rdivide(a.func,b.func); % note the vectorized form
                obj.is_vector=numel(obj.func)>1;
            elseif any(b.func==0)
                % x/0 impossible
                error('dividing by zero not allowed')
            elseif is_one(b)
                % x/1 = x
                obj=a;
            elseif is_zero(a)
                obj.func=0;
            elseif strcmp(a.func,'uminus')
                if strcmp(b.func,'uminus')
                    % (-x) / (-y) = x * y
                    obj=mrdivide(a.args{1},b.args{1});
                else
                    % (-x) / y = -(x*y)
                    obj=-mrdivide(a.args{1},b);
                end
            else
                obj.func='mrdivide';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=mtimes(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            isnum_a=isnumeric(a);
            isnum_b=isnumeric(b);
            if isnum_a && isnum_b
                obj.func=times(a.func,b.func); % note the vectorized form
                obj.is_vector=numel(obj.func)>1;
            elseif strcmp(a.func,'uminus')
                if strcmp(b.func,'uminus')
                    % (-x) * (-y) = x * y
                    obj=mtimes(a.args{1},b.args{1});
                else
                    % (-x) * y = -(x*y)
                    obj=-mtimes(a.args{1},b);
                end
            elseif strcmp(b.func,'uminus')
                % x * (-y) = -(x*y)
                obj=-mtimes(a,b.args{1});
            elseif isnum_a && is_zero(a)
                % 0 * x = 0
                obj=a.prototype;obj.prototype=obj;
                obj.func=0;
                %			obj.is_vector=numel(obj.func)>1;
            elseif isnum_a && is_one(a)
                % 1 * x = x
                obj=b;
            elseif isnum_b && is_zero(b)
                % x * 0 = 0
                obj.func=0;
                %		obj.is_vector=numel(obj.func)>1;
            elseif isnum_b && is_one(b)
                % x * 1 = x
                obj=a;
            else
                obj.func='mtimes';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=ne(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=ne(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='ne';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=normcdf(x,mu,sd)
            if nargin<3
                sd=1;
                if nargin<2
                    mu=0;
                end
            end
            obj=[];
            if strcmp(class(x),'planar'),obj=x.prototype;obj.prototype=obj; end
            if strcmp(class(mu),'planar') && isempty(obj),obj=mu.prototype;obj.prototype=obj; end
            if strcmp(class(sd),'planar') && isempty(obj),obj=sd.prototype;obj.prototype=obj; end
            if ~strcmp(class(x),'planar'),atmp=x; x=obj;x.func=atmp; end
            if ~strcmp(class(mu),'planar'),btmp=mu; mu=obj;mu.func=btmp; end
            if ~strcmp(class(sd),'planar'),ctmp=sd; sd=obj;sd.func=ctmp; end
            if isnumeric(x.func) && isnumeric(mu.func) && isnumeric(sd.func)
                obj.func=normcdf(x.func,mu.func,sd.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='normcdf';
                obj.args={x,mu,sd};
                obj.incidence=x.incidence;
                if ~isempty(mu.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=mu.incidence;
                    else
                        obj.incidence=obj.incidence|mu.incidence;
                    end
                end
                if ~isempty(sd.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=sd.incidence;
                    else
                        obj.incidence=obj.incidence|sd.incidence;
                    end
                end
                obj.is_vector=x.is_vector||mu.is_vector||sd.is_vector;
            end
        end
        function obj=normpdf(x,mu,sd)
            if nargin<3
                sd=1;
                if nargin<2
                    mu=0;
                end
            end
            obj=[];
            if strcmp(class(x),'planar'),obj=x.prototype;obj.prototype=obj; end
            if strcmp(class(mu),'planar') && isempty(obj),obj=mu.prototype;obj.prototype=obj; end
            if strcmp(class(sd),'planar') && isempty(obj),obj=sd.prototype;obj.prototype=obj; end
            if ~strcmp(class(x),'planar'),atmp=x; x=obj;x.func=atmp; end
            if ~strcmp(class(mu),'planar'),btmp=mu; mu=obj;mu.func=btmp; end
            if ~strcmp(class(sd),'planar'),ctmp=sd; sd=obj;sd.func=ctmp; end
            if isnumeric(x.func) && isnumeric(mu.func) && isnumeric(sd.func)
                obj.func=normpdf(x.func,mu.func,sd.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='normpdf';
                obj.args={x,mu,sd};
                obj.incidence=x.incidence;
                if ~isempty(mu.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=mu.incidence;
                    else
                        obj.incidence=obj.incidence|mu.incidence;
                    end
                end
                if ~isempty(sd.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=sd.incidence;
                    else
                        obj.incidence=obj.incidence|sd.incidence;
                    end
                end
                obj.is_vector=x.is_vector||mu.is_vector||sd.is_vector;
            end
        end
        function obj=or(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=or(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='or';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=plus(a,b)
            obj=[];
            if strcmp(class(a),'planar')
                obj=a.prototype;obj.prototype=obj;
            end
            if strcmp(class(b),'planar')
                if isempty(obj)
                    obj=b.prototype;obj.prototype=obj;
                    atmp=a;a=obj; a.func=atmp;
                end
            else
                btmp=b;b=obj; b.func=btmp;
            end
            if isnumeric(a.func) && isnumeric(b.func)
                obj.func=plus(a.func,b.func);
                obj.is_vector=numel(obj.func)>1;
            elseif ischar(b.func) && strcmp(b.func,'uminus')
                % x + (-y) = x - y
                obj=minus(a,b.args{1});
            elseif is_zero(a)
                % 0 + x = x
                obj=b;
            elseif is_zero(b)
                % x + 0 = x
                obj=a;
            elseif isequal(a,b) %% this is going to be expensive ... a.pointer==b.pointer
                % x + x = 2x
                obj=2*a;
            else
                obj.func='plus';
                obj.args={a,b};
                obj.is_vector=a.is_vector||b.is_vector;
                obj.incidence=a.incidence;
                if ~isempty(b.incidence)
                    if isempty(obj.incidence)
                        obj.incidence=b.incidence;
                    else
                        obj.incidence=obj.incidence|b.incidence;
                    end
                end
            end
        end
        function obj=sign(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=sign(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='sign';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=sin(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=sin(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='sin';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=sinh(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=sinh(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='sinh';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=sqrt(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=sqrt(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='sqrt';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=tan(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=tan(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='tan';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=tanh(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=tanh(a.func);
                obj.is_vector=numel(obj.func)>1;
            else
                obj.func='tanh';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=uminus(a)
            obj=a.prototype;obj.prototype=obj;
            if isnumeric(a.func)
                obj.func=uminus(a.func);
                obj.is_vector=numel(obj.func)>1;
            elseif strcmp(a.func,'uminus')
                % Simplify -(-x) in x
                if strcmp(class(a.args{1}),'planar')
                    obj=a.args{1};
                else
                    obj.func=a.args{1};
                end
            else
                obj.func='uminus';
                obj.args={a};
                obj.incidence=a.incidence;
                obj.is_vector=a.is_vector;
            end
        end
        function obj=uplus(a)
            obj=a;
        end
        % bizarre functions
        %------------------
        function obj=power(varargin)
            obj=mpower(varargin{:});
        end
        function obj=rdivide(varargin)
            obj=mrdivide(varargin{:});
        end
        function obj=times(varargin)
            obj=mtimes(varargin{:});
        end
        % utility functions
        %------------------
        varargout=char(varargin)
        varargout=kron(varargin)
        varargout=load_varlist(varargin)
        varargout=set(varargin)
        varargout=get(varargin)
        function flag=isnumeric(a)
            flag=isnumeric(a.func)||islogical(a.func);
        end
        function flag=is_zero(a)
            afunc=a.func;
            flag=(isnumeric(afunc)||islogical(afunc)) && all(afunc==0);
        end
        function flag=is_one(a)
            afunc=a.func;
            flag=(isnumeric(afunc)||islogical(afunc)) && all(afunc==1);
        end
        function obj=intercept_column(obj,pointer)
            if obj.is_vector
                if isnumeric(obj) && numel(obj.func)>1
                    obj.func=obj.func(pointer);
                    obj.is_vector=false;
                elseif ~isempty(obj.args)
                    for iarg=1:numel(obj.args)
                        if ~isa(obj.args{iarg},'planar')
                            continue
                        end
                        obj.args{iarg}=intercept_column(obj.args{iarg},pointer);
                    end
                    % the line below will correct for zeros, ones, etc. as well as
                    % incidences.
                    obj=feval(obj.func,obj.args{:});
                end
            end
        end
        function d=diff(x,wrt,pointer)
            if nargin<3
                pointer=[];
            end
            
            [nrows,ncols]=size(x);
            if nrows>1||ncols>1
                d=cell(nrows,ncols);
                for irow=1:nrows
                    for icol=1:ncols
                        d{irow,icol}=diff(x(irow,icol),wrt);
                    end
                end
                return
            end
            
            if isempty(x.incidence)
                % numbers/vectors and variables which are not part of differentiation
                % automatically receive 0 as derivative
                d=x.prototype;d.prototype=d;d.func=0;
            elseif isempty(x.args)
                % variables that are part of differentiation
                d=x.prototype;d.prototype=d;
                d.func=double(x.incidence(wrt));
                d.is_vector=numel(d.func)>1;
                % why do we need to double it?
                % we need to carry around vectors in order to be able to do
                % higher-order derivatives
            else
                % functions of variables
                if_elseif_flag=strcmp(x.func,'if_elseif');
                if_then_else_flag=strcmp(x.func,'if_then_else');
                nargs=numel(x.args);
                d_args=x.args;
                for iarg=1:nargs
                    if ~isempty(pointer) && isnumeric(x.args{iarg}) && numel(x.args{iarg}.func)>1
                        x.args{iarg}.func=x.args{iarg}.func(pointer);
                        x.args{iarg}.is_vector=false;
                    end
                    if (if_elseif_flag && rem(iarg,2))||(iarg==1 && if_then_else_flag)
                        continue
                    end
                    d_args{iarg}=diff(x.args{iarg},wrt,pointer);
                end
                % compose derivatives
                %--------------------
                switch x.func
                    case 'plus'
                        d=d_args{1}+d_args{2};
                    case 'minus'
                        d=d_args{1}-d_args{2};
                    case 'mtimes'
                        t11 = d_args{1}*x.args{2};
                        t12 = d_args{2}*x.args{1};
                        d=t11+t12;
                    case 'mrdivide'
                        if ~is_zero(d_args{2})
                            t11 = d_args{1}*x.args{2};
                            t12 = d_args{2}*x.args{1};
                            t13 = t11-t12;
                            t14 = x.args{2}^2;
                            d=t13/t14;
                        else
                            d=d_args{1}/x.args{2};
                        end
                    case {'lt','gt','le','ge','eq','ne','or','and','sign'}
                        % derivatives are zero
                        d=x.prototype; d.prototype=d; d.func=0;
                    case 'if_then_else'
                        d=if_then_else(x.args{1},d_args{2},d_args{3});
                    case 'if_elseif'
                        the_args=x.args;
                        the_args(2:2:end)=d_args(2:2:end);
                        d=if_elseif(the_args{:});
                    case 'mpower'
                        t11 = log(x.args{1});
                        t12 = d_args{2}*t11;
                        t13 = d_args{1}*x.args{2};
                        d= t12*x+t13*x.args{1}^(x.args{2}-1);
                        %         t14 = t13/x.args{1}; % <-- creates problems when x.args{1}=0
                        %         t15 = t12+t14;
                        %         d= t15*x.args{1}^x.args{2};
                    case 'max'
                        t11 = x.args{1}>x.args{2};
                        t12 = t11*d_args{1};
                        t13 = 1-t11;
                        t14 = t13*d_args{2};
                        d= t14+t12;
                    case 'min'
                        t11 = x.args{2}>x.args{1};
                        t12 = t11*d_args{1};
                        t13 = 1-t11;
                        t14 = t13*d_args{2};
                        d= t14+t12;
                    case 'uminus'
                        d= -d_args{1};
                    case 'uplus'
                        d= d_args{1};
                    case 'exp'
                        d= d_args{1}*x;
                    case 'log'
                        d= d_args{1}/x.args{1};
                    case 'log10'
                        t11 = exp(1);
                        t12 = log10(t11);
                        t13 = d_args{1}/x.args{1};
                        d= t12*t13;
                    case 'cos'
                        t11 = sin(x.args{1});
                        t12 = -t11;
                        d= d_args{1}*t12;
                    case 'sin'
                        t11 = cos(x.args{1});
                        d= d_args{1}*t11;
                    case 'tan'
                        d= d_args{1}*(1+x^2);
                    case 'cot'
                        d= -d_args{1}*(1-x^2); % <--- d=-d_args{1}/sin(x.args{1})^2;
                    case 'acos'
                        d= -d_args{1}/sqrt(1-x^2);
                    case 'asin'
                        t11 = cos(x);
                        d= d_args{1}/t11;
                    case 'atan'
                        t11 = x.args{1}^2;
                        t12 = 1+t11;
                        d= d_args{1}/t12;
                    case 'cosh'
                        t11 = sinh(x.args{1});
                        d= d_args{1}*t11;
                    case 'sinh'
                        t11 = cosh(x.args{1});
                        d= d_args{1}*t11;
                    case 'tanh'
                        d= d_args{1}*(1-x^2);
                    case 'acosh'
                        t11 = sinh(x);
                        d= d_args{1}/t11;
                    case 'asinh'
                        t11 = cosh(x);
                        d= d_args{1}/t11;
                    case 'atanh'
                        t11 = x.args{1}^2;
                        t12 = 1-t11;
                        d= d_args{1}*t12;
                    case 'sqrt'
                        t11 = 2*x;
                        d= d_args{1}/t11;
                    case 'abs'
                        t11 = sign(x.args{1});
                        d= t11*d_args{1};
                    case 'erf'
                        % x^2
                        t11 = mpower(x.args{1},2);
                        % exp(x^2)
                        t12 =  exp(t11);
                        % sqrt(pi)
                        t11 = sqrt(pi);
                        % sqrt(pi)*exp(x^2)
                        t13 = t11*t12;
                        % 2/(sqrt(pi)*exp(x^2));
                        t14 = 2/t13;
                        % (2/(sqrt(pi)*exp(x^2)))*dx;
                        d= t14*d_args{1};
                    case 'normcdf'
                        d=(d_args{1}+d_args{2}+d_args{3})*normpdf(x.args{1},x.args{2},x.args{3});
                    case 'normpdf'
                        y=(x.args{1}-x.args{2})/x.args{3};
                        ss= d_args{3}/x.args{3};
                        d=((ss*y-(d_args{1}-d_args{2})/x.args{3})*y-ss)*x;
                end
            end
        end
    end
    methods(Static)
        function deriv=derivatives2functions(deriv,args,compact,optimize)
            if nargin<4
                optimize=false;
                if nargin<3
                    compact=false;
                end
            end
            if ischar(args)
                args=cellstr(args);
            end
            args=args(:)';
            % put into analytical form
            %--------------------------
            deriv.derivatives=parser.analytical_symbolic_form(deriv.derivatives,args,'analytic');
            % remove unnecessary parentheses
            %-------------------------------
            if optimize
                word='\w+';
                word_par='\w+\(\d+\)';
                deriv.derivatives=regexprep(deriv.derivatives,...
                    ['(?<!\w+)(\()(',word,'|',word_par,')(\))'],'$2');
            end
            % build the functions
            %--------------------
            args=cell2mat(strcat(args,','));
            main_string=['@(',args(1:end-1),')'];
            nrows=deriv.size(1);
            maxcols=deriv.maxcols;
            for irow=1:nrows
                if compact
                    dd='';
                    dindex=[];
                end
                for icol=1:maxcols
                    if isempty(deriv.derivatives{irow,icol})
                        break
                    end
                    icol_deriv=deriv.derivatives{irow,icol};
                    if compact
                        dd=[dd,icol_deriv,','];
                        dindex=[dindex,deriv.map{irow,icol}];
                    else
                        deriv.derivatives{irow,icol}=str2func([main_string,icol_deriv]);
                    end
                end
                if compact && ~isempty(dd)
                    deriv.derivatives{irow,1}=str2func([main_string,'[',dd(1:end-1),']']);
                    deriv.map{irow,1}=dindex;
                end
            end
            if compact && ~isempty(deriv.derivatives)
                deriv.derivatives=deriv.derivatives(:,1);
                deriv.map=deriv.map(:,1);
                deriv.maxcols=1;
            end
        end
        function c=print(deriv,long)
            if nargin<2
                long=false;
            end
            nrows=deriv.size(1);
            c=cell(nrows,deriv.maxcols);
            cmap=c;
            for irow=1:nrows
                this=deriv.derivatives{irow};
                ncols=numel(this);
                for icol=1:ncols
                    if isempty(this(icol).location)
                        continue
                    end
                    c{irow,icol}=char(this(icol),long);
                    cmap{irow,icol}=this(icol).location;
                end
            end
            deriv.derivatives=c;
            c=deriv;
            c.map=cmap;
        end
        function var_list=initialize(var_list,wrt_list)
            if ischar(var_list)
                var_list=cellstr(var_list);
            end
            if ischar(wrt_list)
                wrt_list=cellstr(wrt_list);
            end
            nwrt=numel(wrt_list);
            nvlist=numel(var_list);
            if numel(unique(wrt_list))~=nwrt
                error('repeated variable names in wrt_list')
            end
            if numel(unique(var_list))~=nvlist
                error('repeated variable names in var_list')
            end
            proto_incidence=false(1,nwrt);
            proto_=planar();
            for ivar=1:nvlist
                vname=var_list{ivar};
                loc=find(strcmp(vname,wrt_list));
                var_list{ivar}=planar(vname);
                if ~isempty(loc)
                    incid_=proto_incidence;
                    incid_(loc)=true;
                    %         var_list{ivar}.incid_=sparse(incid_);
                    var_list{ivar}=set(var_list{ivar},'incidence',sparse(incid_));
                end
                var_list{ivar}=set(var_list{ivar},'prototype',proto_);
            end
        end
        function derivs=differentiate(eqtns,nwrt,order,partitions,verbose,early_intercept)
            if nargin<6
                early_intercept=true;
                if nargin<5
                    verbose=false;
                    if nargin<4
                        partitions=[];
                    end
                end
            end
            wrt_wise=true;
            
            neqtns=size(eqtns,1);
            
            [labels,low,high,nlabels]=set_partitions();

            init_vector=repmat({planar.empty(0)},neqtns,1);
            derivs=struct('size',{},'derivatives',{},'maxcols',{},...
                'nnz_derivs',{},'map',{},'partitions',{});
            combos=[];
            old_ncols=1;
            MainGrid=[];
            LabelGrid=[];
            for oo=1:order
                if verbose
                    tic
                end
                % build the grids
                %----------------
                MainGrid=utils.gridfuncs.build_grid(MainGrid,nwrt);
                LabelGrid=utils.gridfuncs.build_grid(LabelGrid,nlabels);
                
                % prepare permutations
                %---------------------
                proto_permutation=cell2mat(utils.gridfuncs.mypermutation(1:oo));

                % compute all combinations and hash them up
                %------------------------------------------
                [combos,iG_DerivLocs,ncols]=store_combinations(combos);
                
                % store the partitions
                %---------------------
                DerivPartitions=store_derivative_partitions();
                
                % differentiate
                %--------------
                d=init_vector;
                maxcols=0;
                nnz_derivs=0;
                for ieqtn=1:neqtns % number of rows
                    maxcols_ielt=0;
                    for ielt=1:numel(eqtns{ieqtn}) % number of elements in the row
                        oldwrt=eqtns{ieqtn}(ielt).lineage;
                        if ~isempty(oldwrt)
                            oldwrt=oldwrt{end};
                        end
                        if ~early_intercept
                            obj=eqtns{ieqtn}(ielt);
                            candidates=find(obj.incidence);
                        end
                        for icol=1:max(1,numel(oldwrt)) % number of derivatives taken in earlier round
                            if early_intercept
                                obj=intercept_column(eqtns{ieqtn}(ielt),icol);
                                if ~isempty(eqtns{ieqtn}(ielt).lineage) && ~isempty(eqtns{ieqtn}(ielt).lineage{end})
                                    obj.lineage=[eqtns{ieqtn}(ielt).lineage(1:end-1),{eqtns{ieqtn}(ielt).lineage{end}(icol)}];
                                end
                                candidates=find(obj.incidence);
                                pointer=[];
                            else
                                pointer=icol;
                            end
                            wrt=candidates;
                            if ~isempty(oldwrt) % <--- oo>1
                                wrt(wrt<oldwrt(icol))=[];
                                lineage__=[obj.lineage(1:end-1),{oldwrt(icol),wrt}];
                            else
                                lineage__={wrt};
                            end
                            d0=diff(obj,wrt,pointer);
                            d0.lineage=lineage__;
                            nderivs=numel(d0.lineage{end});
                            % locate each derivative in the reduced matrix
                            %---------------------------------------------
                            former=cell2mat(d0.lineage(1:end-1));
                            index_=nan(1,nderivs);
                            for id=1:nderivs
                                index_(id)=iG_DerivLocs(...
                                    utils.gridfuncs.locate_permutation([former,lineage__{end}(id)],nwrt,wrt_wise)...
                                    );
                            end
                            d0.location=index_;
                            % store the derivative if not zero
                            %---------------------------------
                            if ~isempty(index_)
                                maxcols_ielt=maxcols_ielt+1;
                                d{ieqtn}=[d{ieqtn},d0];
                                nnz_derivs=nnz_derivs+numel(index_);
                            end
                        end
                        maxcols=max(maxcols,maxcols_ielt);
                    end
                end
                derivs(oo)=struct('size',{[neqtns,ncols]},'derivatives',{d},...
                    'maxcols',maxcols,'nnz_derivs',nnz_derivs,'map',[],...
                    'partitions',DerivPartitions);

                % spit out the time it took to compute
                %-------------------------------------
                if verbose
                    fprintf(1,'differentiation at order %0.0f done in %0.4f seconds\n',oo,toc);
                end
                
                % next round
                %-----------
                eqtns=d;
            end
            function [newcombos,iG_DerivLocs,iter]=store_combinations(oldcombos)
                iG_DerivLocs=nan(1,nwrt^oo);
                iter=0;
                newcombos=cell(1,nwrt*old_ncols);
                for icol_=1:old_ncols
                    if isempty(oldcombos)
                        main=[];
                        start=1;
                    else
                        main=oldcombos{icol_};
                        start=main(end);
                    end
                    for irt=start:nwrt
                        iter=iter+1;
                        index=[main,irt];
                        newcombos{iter}=index;
                        % find the location in iG
                        %------------------------
                        pp=index(proto_permutation);
                        extended_locs=utils.gridfuncs.locate_permutation(pp,nwrt,wrt_wise);
                        % remove redundancies
                        extended_locs=unique(extended_locs);
                        iG_DerivLocs(extended_locs)=iter;
                   end
                end
                old_ncols=iter;
                newcombos=newcombos(1:iter);
            end
            function  DerivPartitions=store_derivative_partitions()
                iG=MainGrid(:,end:-1:1); % flip the grid to put it into the wrt form
                test=true;
                if test
                    iG_DerivLocsShort=iG_DerivLocs;
                    position=true(nwrt^oo,1);
                    DerivPartitions=struct();
                    for irow=1:size(LabelGrid,1)
                        for icol_=1:oo
                            position=position & iG(:,icol_)>low(LabelGrid(irow,icol_)) & iG(:,icol_)<high(LabelGrid(irow,icol_));
                        end
                        xxx=cell2mat(labels(LabelGrid(irow,:)));
                        DerivPartitions.(xxx)=iG_DerivLocsShort(position);
                        iG=iG(~position,:); % <---iG(position,:)=[];
                        iG_DerivLocsShort=iG_DerivLocsShort(~position);
                        position=true(size(iG,1),1);
                    end
                else
                    id0=true(nwrt^oo,1);
                    for irow=1:size(LabelGrid,1)
                        position=id0;
                        for icol_=1:oo
                            position=position & iG(:,icol_)>low(LabelGrid(irow,icol_)) & iG(:,icol_)<high(LabelGrid(irow,icol_));
                        end
                        xxx=cell2mat(labels(LabelGrid(irow,:)));
                        DerivPartitions{oo}.(xxx)=iG_DerivLocs(position);
                    end
                end
            end
            function [labels,low,high,nparts]=set_partitions()
                if isempty(partitions)
                    partitions={'d';nwrt};
                end
                parts=cell2mat(partitions(2,:));
                if sum(parts)~=nwrt
                    error('the sum of the partitions should be equal to the number of variables to be differentiated')
                end
                labels=partitions(1,:);
                nparts=numel(parts);
                low=nan(1,nparts);
                high=nan(1,nparts);
                tt=cumsum(parts);
                tt=[0,tt];
                for ipart=1:nparts
                    low(ipart)=tt(ipart);
                    high(ipart)=tt(ipart+1)+1;
                end
            end
        end
    end
end
