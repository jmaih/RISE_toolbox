classdef rise_sad < handle
    properties
        name
        args
        ref
        ncalls=0
    end
    methods
        function obj=rise_sad(name,arg)
            if nargin
                if isa(name,'rise_sad')
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
                            arg{iarg}=rise_sad(arg{iarg});
                        end
                    end
                    obj.args=arg;
                end
            end
        end
        function obj=plus(u,v),obj=rise_sad('plus',{u,v}); end
        function obj=uplus(u),obj=rise_sad('uplus',{u}); end
        function obj=minus(u,v),obj=rise_sad('minus',{u,v}); end
        function obj=uminus(u),obj=rise_sad('uminus',{u}); end
        function obj=times(u,v),obj=rise_sad('times',{u,v}); end
        function obj=mtimes(u,v),obj=rise_sad('mtimes',{u,v}); end
        function obj=power(u,v),obj=rise_sad('power',{u,v}); end
        function obj=mpower(u,v),obj=rise_sad('mpower',{u,v}); end
        function obj=rdivide(u,v),obj=rise_sad('rdivide',{u,v}); end
        function obj=mrdivide(u,v),obj=rise_sad('mrdivide',{u,v}); end
        function obj=ldivide(u,v),obj=rise_sad('ldivide',{u,v}); end
        function obj=mldivide(u,v),obj=rise_sad('mldivide',{u,v}); end
        function obj=exp(u),obj=rise_sad('exp',{u}); end
        function obj=log(u),obj=rise_sad('log',{u}); end
        function obj=log10(u),obj=rise_sad('log10',{u}); end
        function obj=cos(u),obj=rise_sad('cos',{u}); end
        function obj=acos(u),obj=rise_sad('acos',{u}); end
        function obj=cosh(u),obj=rise_sad('cosh',{u}); end
        function obj=sin(u),obj=rise_sad('sin',{u}); end
        function obj=asin(u),obj=rise_sad('asin',{u}); end
        function obj=sinh(u),obj=rise_sad('sinh',{u}); end
        function obj=tan(u),obj=rise_sad('tan',{u}); end
        function obj=atan(u),obj=rise_sad('atan',{u}); end
        function obj=tanh(u),obj=rise_sad('tanh',{u}); end
        function obj=min(u,v),obj=rise_sad('min',{u,v}); end
        function obj=max(u,v),obj=rise_sad('max',{u,v}); end
        function obj=sum(u,v),obj=rise_sad('sum',{u,v}); end
        function obj=normpdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=rise_sad('normpdf',{u,v,w});
        end
        function obj=normcdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=rise_sad('normcdf',{u,v,w});
        end
        function obj=abs(u),obj=rise_sad('abs',{u}); end
        function obj=isreal(u),obj=rise_sad('isreal',{u}); end
        function obj=sqrt(u),obj=rise_sad('sqrt',{u}); end
        function obj=norm(u),obj=rise_sad('norm',{u}); end
        function obj=gt(u,v),obj=rise_sad('gt',{u,v}); end
        function obj=ge(u,v),obj=rise_sad('ge',{u,v}); end
        function obj=lt(u,v),obj=rise_sad('lt',{u,v}); end
        function obj=le(u,v),obj=rise_sad('le',{u,v}); end
        function obj=sign(u),obj=rise_sad('sign',{u}); end
        function string=char(obj,unravel,isparent)
            if nargin<3
                isparent=false;
                if nargin<2
                    unravel=false;
                end
            end
            n=numel(obj);
            if n>1
                string=cell(size(obj));
                for i0=1:n
                    string{i0}=char(obj(i0),unravel,isparent);
                end
                return
            end
            % char itself is already taken care of
            args_=reprocess_arguments(obj.args);
            if isa(obj.name,'double')
                string=num2str(obj.name,10);
            elseif isempty(args_) % variable
                string=obj.name; % this should be a char
            elseif ~unravel && ~isempty(obj.ref) && ~isparent
                string=obj.ref;
            else
                switch obj.name
                    case {'plus','minus','times','power'}
                        string=rise_cat(mychar(args_{1},unravel),mychar(args_{2},unravel),obj.name);
                    case {'uplus','uminus'}
                        string=rise_cat('0',mychar(args_{1},unravel),obj.name(2:end));
                    case {'mtimes','mpower','rdivide'}
                        string=rise_cat(mychar(args_{1},unravel),mychar(args_{2},unravel),obj.name(2:end));
                    case {'mrdivide'}
                        string=rise_cat(mychar(args_{1},unravel),mychar(args_{2},unravel),obj.name(3:end));
                    case {'min','max','gt','lt','ge','le'}
                        string=[obj.name,'(',mychar(args_{1},unravel),',',mychar(args_{2},unravel),')'];
                    case {'ldivide','mldivide'}
                        string=rise_cat(mychar(args_{2},unravel),mychar(args_{1},unravel),'divide');
                    case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                            'tan','atan','tanh','abs','sqrt','isreal','sign'}
                        string=[obj.name,'(',mychar(args_{1},unravel),')'];
                    case {'normpdf','normcdf'}
                        string=[obj.name,'(',mychar(args_{1},unravel),',',mychar(args_{2},unravel),',',mychar(args_{3},unravel),')'];
                end
            end
            if isparent
                string=[obj.ref,'=',string,';'];
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
        function flag=eq(obj1,obj2)
            siz=size(obj1);
            siz2=size(obj2);
            if isequal(siz,siz2)
            elseif max(siz)==1
                obj1=repmat(obj1,siz2);
                siz=siz2;
            elseif max(siz2)==1
                obj2=repmat(obj2,siz);
            else
                error('wrong sizes of inputs')
            end
            flag=false(siz);
            for ii=1:prod(siz)
                if isequal(obj1(ii).name,obj2(ii).name)
                    nargs=numel(obj1(ii).args);
                    flag(ii)=true;
                    if nargs
                        for iarg=1:nargs
                            eq_args=eq(obj1(ii).args{iarg},obj2(ii).args{iarg});
                            if ~eq_args
                                flag(ii)=false;
                                break
                            end
                        end
                    end
                end
            end
        end
        function [this,references,tank]=diff(obj,wrt,references,tank)
            % check that I do not need to output the object thanks to the
            % handle property
            if nargin<4
                tank=[];
                if nargin<3
                    references=[];
                end
            end
            tank=consolidate(obj,tank);
            this=mydiff(obj,wrt);
%             tank=consolidate(this,tank);
%             tank=update_flags(obj,tank);
%             tank=update_flags(this,tank);
            references=collect_references(obj,references);
%             references=collect_references(this,references);
        end
    end
    methods(Access=private)
        function tank=update_flags(tree,tank)
            update_flags_intern(tree,1);
            update_flags_intern(tree,2);
            function update_flags_intern(tree,pass)
                for it=1:numel(tree)
                    if ~isempty(tree(it).args)
                        for iarg=1:numel(tree(it).args)
                            update_flags_intern(tree(it).args{iarg},pass);
                        end
                        if ~isempty(tree(it).ref)
                            location=str2double(tree(it).ref(3:end));
                            calls=tank{3}(location);
                            if pass==1
                                tank{3}(location)=calls+tree(it).ncalls;
                            elseif pass==2
                                tree(it).ncalls=calls;
%                                 if calls<2
%                                     tree(it).ref=[];
%                                 end
                            else
                                error('pass must be either 1 or 2')
                            end
                        end
                    end
                end
            end
        end
        function tank=consolidate(tree,tank)
            % vectorizes and consolidates a tree vector... 
            % I should check that I do not need to call the output because
            % of the handle property
            if nargin<2
                tank=[];
            end
            if isempty(tank)
                tank=containers.Map();
            end
            for it=1:numel(tree)
                if ~isempty(tree(it).args)
                    string=char(tree(it),true);
                    if rise_isa(string,'atom')
                        tree(it)=rise_sad(string);
                    elseif isKey(tank,string)
                        tree(it)=tank(string);
                    else
                        for iarg=1:numel(tree(it).args)
                            tank=consolidate(tree(it).args{iarg},tank);
                        end
                        n=tank.Count+1;
                        tree(it).ref=['T_',int2str(n)];
                        tank(string)=tree(it);
                    end
                end
            end
        end
        function update_reference(varargin)
%             not a variable and not a constant
            for iarg=1:length(varargin)
                if ~isempty(varargin{iarg}.ref) && ~isempty(varargin{iarg}.args)
                    varargin{iarg}.ncalls=varargin{iarg}.ncalls+1;
                end
            end
        end
        function references=collect_references(tree,references)
            if nargin<2
                references=[];
            end
            unravel=false;
            for it=1:numel(tree)
                if ~isempty(tree(it).args)
                    isparent=tree(it).ncalls>=2;% <---is_atom(char(tree(it),false));
%                     isparent=~isempty(tree(it).ref) && strcmp(tree(it).ref(1),'T');% <---is_atom(char(tree(it),false));
                    for iarg=1:numel(tree(it).args)
                        references=collect_references(tree(it).args{iarg},references);
                    end
                    if isparent
                        new_item=char(tree(it),unravel,isparent);
                        if ~any(strcmp(new_item,references))
                            references=[references;{new_item}]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        
        function this=mydiff(obj,wrt)
            this=rise_sad.empty(0,0);
            for iobj=1:numel(obj)
                if isnumeric(obj(iobj).name)
                    this(iobj,1)=rise_sad(0);
                elseif isempty(obj(iobj).args)
                    if isequal(obj(iobj).name,wrt.name)
                        this(iobj,1)=rise_sad(1);
                    else
                        this(iobj,1)=rise_sad(0);
                    end
                else
                    u=rise_sad(obj(iobj).args{1});
                    n_args=numel(obj(iobj).args);
                    if n_args>1
                        v=rise_sad(obj(iobj).args{2});
                        if n_args>2
                            w=rise_sad(obj(iobj).args{3});
                        end
                    end
                    if iscell(wrt)
                        for ic=1:numel(wrt)
                            this(iobj,ic)=differentiation_engine(iobj,wrt{ic});
                        end
                    elseif isa(wrt,'rise_sad')
                        this(iobj,1)=differentiation_engine(iobj,wrt);
                    else
                        error([mfilename,':: second argument must be rise_sad object'])
                    end
                end
            end
            function this=differentiation_engine(index,wrt)
                du=mydiff(u,wrt);
                switch obj(index).name
                    case {'gt','ge','lt','le','sign'}
                        this=rise_sad(0);
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
                        update_reference(u,v);
                    case {'power','mpower'}
                        this=v*du*u^(v-1);
                        if ~isnumeric(v.name)
                            dv=mydiff(v,wrt);
                            this=this+dv*log(u)*u^v;
                        end
                        update_reference(u,v);
                    case {'rdivide','mrdivide'}
                        dv=mydiff(v,wrt);
                        this=(du*v-dv*u)/v^2; % only when v~=0
                        update_reference(u,v);
                    case {'ldivide','mldivide'}
                        dv=mydiff(v,wrt);
                        this=(u*dv-v*du)/u^2; % only when u~=0
                        update_reference(u,v);
                    case 'exp'
                        this=du*obj(index);
                        update_reference(obj(index));
                    case 'log'
                        this=du/u;
                        update_reference(u);
                    case 'log10'
                        this=(du/u)/log(10);
                        update_reference(u);
                    case 'cos'
                        this=-du*sin(u);
                        update_reference(u);
                    case 'acos'
                        this=-du/sqrt(1-u^2);
                        update_reference(u);
                    case 'cosh'
                        this=du*sinh(u);
                        update_reference(u);
                    case 'sin'
                        this=du*cos(u);
                        update_reference(u);
                    case 'asin'
                        this=du/sqrt(1-u^2);
                        update_reference(u);
                    case 'sinh'
                        this=du*cosh(u);
                        update_reference(u);
                    case 'tan'
                        this=du/(cos(u))^2;
                        update_reference(u);
                    case 'atan'
                        this=du/(1+u^2);
                        update_reference(u);
                    case 'tanh'
                        this=du/(1-u^2);
                        update_reference(u);
                    case 'min'
                        muv=u<v;
                        this=muv*u+(1-muv)*v;
                        update_reference(u,v);
                    case 'max'
                        muv=u>v;
                        this=muv*u+(1-muv)*v;
                        update_reference(u,v);
                    case 'sum'
                        this=sum(du);
                        update_reference(u);
                    case 'normpdf'
                        this=-du/w*(u-v)/w*obj(index);
                        updtate_reference(u,v,w);
                    case 'normcdf'
                        this=du*normpdf(u,v,w);
                        update_reference(u,v,w);
                    case 'abs'
                        this=du*(-u<0+u>0);
                        update_reference(u);
                    case 'isreal'
                        this=isreal(u)*du;
                        update_reference(u);
                    case 'sqrt'
                        this=du/(2*sqrt(obj(index)));
                        update_reference(obj(index));
                    case 'norm' % this would not work!
                        this=sum(u.*du)/norm(u);
                        update_reference(u);
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

function cc=mychar(obj,unravel)
if ischar(obj)
    cc=obj;
else
    cc=char(obj,unravel);
end
end

