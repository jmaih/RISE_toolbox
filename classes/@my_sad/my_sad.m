classdef my_sad < handle
    properties
        name
        args
    end
    methods
        function obj=my_sad(name,arg)
            if nargin
                if isa(name,'my_sad')
                    obj=name;
                else
                    obj.name=name;
                    if nargin<2
                        arg=[];
                    end
                    if isempty(arg)
                        arg={};
                    end
                    if ~iscell(arg)
                        arg={arg};
                    end
                    for iarg=1:numel(arg)
                        if isnumeric(arg{iarg})
                            arg{iarg}=my_sad(arg{iarg});
                        end
                    end
                    obj.args=arg;
                end
%             else
%                 mymethods=methods(mfilename);
%                 disp(mymethods)
            end
        end
        function obj=plus(u,v),obj=my_sad('plus',{u,v}); end
        function obj=uplus(u),obj=my_sad('uplus',{u}); end
        function obj=minus(u,v),obj=my_sad('minus',{u,v}); end
        function obj=uminus(u),obj=my_sad('uminus',{u}); end
        function obj=times(u,v),obj=my_sad('times',{u,v}); end
        function obj=mtimes(u,v),obj=my_sad('mtimes',{u,v}); end
        function obj=power(u,v),obj=my_sad('power',{u,v}); end
        function obj=mpower(u,v),obj=my_sad('mpower',{u,v}); end
        function obj=rdivide(u,v),obj=my_sad('rdivide',{u,v}); end
        function obj=mrdivide(u,v),obj=my_sad('mrdivide',{u,v}); end
        function obj=ldivide(u,v),obj=my_sad('ldivide',{u,v}); end
        function obj=mldivide(u,v),obj=my_sad('mldivide',{u,v}); end
        function obj=exp(u),obj=my_sad('exp',{u}); end
        function obj=log(u),obj=my_sad('log',{u}); end
        function obj=log10(u),obj=my_sad('log10',{u}); end
        function obj=cos(u),obj=my_sad('cos',{u}); end
        function obj=acos(u),obj=my_sad('acos',{u}); end
        function obj=cosh(u),obj=my_sad('cosh',{u}); end
        function obj=sin(u),obj=my_sad('sin',{u}); end
        function obj=asin(u),obj=my_sad('asin',{u}); end
        function obj=sinh(u),obj=my_sad('sinh',{u}); end
        function obj=tan(u),obj=my_sad('tan',{u}); end
        function obj=atan(u),obj=my_sad('atan',{u}); end
        function obj=tanh(u),obj=my_sad('tanh',{u}); end
        function obj=min(u,v),obj=my_sad('min',{u,v}); end
        function obj=max(u,v),obj=my_sad('max',{u,v}); end
        function obj=sum(u,v),obj=my_sad('sum',{u,v}); end
        function obj=normpdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=my_sad('normpdf',{u,v,w});
        end
        function obj=normcdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=my_sad('normcdf',{u,v,w});
        end
        function obj=abs(u),obj=my_sad('abs',{u}); end
        function obj=isreal(u),obj=my_sad('isreal',{u}); end
        function obj=sqrt(u),obj=my_sad('sqrt',{u}); end
        function obj=norm(u),obj=my_sad('norm',{u}); end
        function obj=gt(u,v),obj=my_sad('gt',{u,v}); end
        function obj=ge(u,v),obj=my_sad('ge',{u,v}); end
        function obj=lt(u,v),obj=my_sad('lt',{u,v}); end
        function obj=le(u,v),obj=my_sad('le',{u,v}); end
        function obj=sign(u),obj=my_sad('sign',{u}); end
        function string=char(obj)
            n=numel(obj);
            if n>1
                string=cell(size(obj));
                for i0=1:n
                    string{i0}=char(obj(i0));
                end
                return
            end
            % char itself is already taken care of
            args_=reprocess_arguments(obj.args);
            if isa(obj.name,'double')
                string=num2str(obj.name,10);
            elseif isempty(args_) % variable
                string=obj.name; % this should be a char
            else
                switch obj.name
                    case {'plus','minus','times','power'}
                        string=rise_cat(char(args_{1}),char(args_{2}),obj.name);
                    case {'uplus','uminus'}
                        string=rise_cat('0',char(args_{1}),obj.name(2:end));
                    case {'mtimes','mpower','rdivide'}
                        string=rise_cat(char(args_{1}),char(args_{2}),obj.name(2:end));
                    case {'mrdivide'}
                        string=rise_cat(char(args_{1}),char(args_{2}),obj.name(3:end));
                    case {'min','max','gt','lt','ge','le'}
                        string=[obj.name,'(',char(args_{1}),',',char(args_{2}),')'];
                    case {'ldivide','mldivide'}
                        string=rise_cat(char(args_{2}),char(args_{1}),'divide');
                    case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                            'tan','atan','tanh','abs','sqrt','isreal','sign'}
                        string=[obj.name,'(',char(args_{1}),')'];
                    case {'normpdf','normcdf'}
                        string=[obj.name,'(',char(args_{1}),',',char(args_{2}),',',char(args_{3}),')'];
                    otherwise
                         string=obj.name;
                         for iarg=1:numel(args_)
                             if iarg==1
                                 string=[string,'(']; %#ok<*AGROW>
                             end
                             string=[string,',',char(args_{iarg})];
                             if iarg==numel(args_)
                                 string=[string,')'];
                             end
                         end
                end
            end
        end
        function print(tree)
            for it=1:numel(tree)
                %                 if ~isempty(tree(it).ref)
                %                     fprintf(1,'%s\n',[tree(it).ref,' ---> ',tree(it).name]);
                fprintf(1,'%s\n',[' ---> ',tree(it).name]);
                %                 end
                tree(it).args=reprocess_arguments(tree(it).args);
                for iarg=1:numel(tree(it).args)
                    if ~ischar(tree(it).args{iarg})
                        print(tree(it).args{iarg})
                    end
                end
            end
        end
        function this=diff(obj,wrt)
            this=mydiff(obj,wrt);
        end
    end
    methods(Access=private)        
        function this=mydiff(obj,wrt)
            this=my_sad.empty(0,0);
            for iobj=1:numel(obj)
                if isnumeric(obj(iobj).name)
                    this(iobj,1)=my_sad(0);
                elseif isempty(obj(iobj).args)
                    if isequal(obj(iobj).name,wrt.name)
                        this(iobj,1)=my_sad(1);
                    else
                        this(iobj,1)=my_sad(0);
                    end
                else
                    u=my_sad(obj(iobj).args{1});
                    n_args=numel(obj(iobj).args);
                    if n_args>1
                        v=my_sad(obj(iobj).args{2});
                        if n_args>2
                            w=my_sad(obj(iobj).args{3});
                        end
                    end
                    if iscell(wrt)
                        for ic=1:numel(wrt)
                            this(iobj,ic)=differentiation_engine(iobj,wrt{ic});
                        end
                    elseif isa(wrt,'my_sad')
                        this(iobj,1)=differentiation_engine(iobj,wrt);
                    elseif isa(wrt,'char')
                        this(iobj,1)=differentiation_engine(iobj,my_sad(wrt));
                    else
                        error([mfilename,':: second argument must be my_sad object'])
                    end
                end
            end
            function this=differentiation_engine(index,wrt)
                du=mydiff(u,wrt);
                switch obj(index).name
                    case {'gt','ge','lt','le','sign'}
                        this=my_sad(0);
                    case 'plus'
                        dv=mydiff(v,wrt);
                        this=du+dv;
                    case 'uplus'
                        this=du;
                    case 'minus'
                        dv=mydiff(v,wrt);
                        this=du-dv;
                    case 'uminus'
                        this=-du;
                    case {'mtimes','times'}
                        dv=mydiff(v,wrt);
                        upv=du*v;
                        vpu=dv*u;
                        this=upv+vpu;
                    case {'power','mpower'}
                        this=v*du*u^(v-1);
                        if ~isnumeric(v.name)
                            dv=mydiff(v,wrt);
                            this=this+dv*log(u)*u^v;
                        end
                    case {'rdivide','mrdivide'}
                        dv=mydiff(v,wrt);
                        this=(du*v-dv*u)/v^2; % only when v~=0
                    case {'ldivide','mldivide'}
                        dv=mydiff(v,wrt);
                        this=(u*dv-v*du)/u^2; % only when u~=0
                    case 'exp'
                        this=du*obj(index);
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
                        update_reference(u);
                    case 'min'
                        muv=u<v;
                        this=muv*u+(1-muv)*v;
                    case 'max'
                        muv=u>v;
                        this=muv*u+(1-muv)*v;
                    case 'sum'
                        this=sum(du);
                    case 'normpdf'
                        this=-du/w*(u-v)/w*obj(index);
                    case 'normcdf'
                        this=du*normpdf(u,v,w);
                        update_reference(u,v,w);
                    case 'abs'
                        this=du*(-u<0+u>0);
                    case 'isreal'
                        this=isreal(u)*du;
                    case 'sqrt'
                        this=du/(2*sqrt(obj(index)));
                    case 'norm' % this would not work!
                        this=sum(u.*du)/norm(u);
                    otherwise
                        error([obj(index).name,' is unknown type of operator'])
                end
            end
        end
    end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
    end
end


function args=reprocess_arguments(args)
for ii=1:numel(args)
    if isnumeric(args{ii})
        args{ii}=num2str(args{ii},10);
    end
end
end


