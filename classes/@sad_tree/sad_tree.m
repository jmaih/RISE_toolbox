classdef sad_tree < handle
    properties
        name
        args
        ref
        ncalls=0
    end
    methods
        function obj=sad_tree(name,arg,initialize)
            if nargin
                if isa(name,'sad_tree')
                    obj=name;
                else
                    if nargin<3
                        initialize=[];
                        if nargin<2
                            arg=[];
                        end
                    end
                    if isempty(initialize),initialize=false;end
                    if isempty(arg),arg={}; end
                    obj.name=name;
                    if ~iscell(arg)
                        arg={arg};
                    end
                    obj.args=arg;
                    if initialize
                        setappdata(0,'sad_tree_formulae',{})
                    end
                end
            end
        end
        function obj=plus(u,v),obj=sad_tree('plus',{u,v}); end
        function obj=uplus(u),obj=sad_tree('uplus',{u}); end
        function obj=minus(u,v),obj=sad_tree('minus',{u,v}); end
        function obj=uminus(u),obj=sad_tree('uminus',{u}); end
        function obj=times(u,v),obj=sad_tree('times',{u,v}); end
        function obj=mtimes(u,v),obj=sad_tree('mtimes',{u,v}); end
        function obj=power(u,v),obj=sad_tree('power',{u,v}); end
        function obj=mpower(u,v),obj=sad_tree('mpower',{u,v}); end
        function obj=rdivide(u,v),obj=sad_tree('rdivide',{u,v}); end
        function obj=mrdivide(u,v),obj=sad_tree('mrdivide',{u,v}); end
        function obj=ldivide(u,v),obj=sad_tree('ldivide',{u,v}); end
        function obj=mldivide(u,v),obj=sad_tree('mldivide',{u,v}); end
        function obj=exp(u),obj=sad_tree('exp',{u}); end
        function obj=log(u),obj=sad_tree('log',{u}); end
        function obj=log10(u),obj=sad_tree('log10',{u}); end
        function obj=cos(u),obj=sad_tree('cos',{u}); end
        function obj=acos(u),obj=sad_tree('acos',{u}); end
        function obj=cosh(u),obj=sad_tree('cosh',{u}); end
        function obj=sin(u),obj=sad_tree('sin',{u}); end
        function obj=asin(u),obj=sad_tree('asin',{u}); end
        function obj=sinh(u),obj=sad_tree('sinh',{u}); end
        function obj=tan(u),obj=sad_tree('tan',{u}); end
        function obj=atan(u),obj=sad_tree('atan',{u}); end
        function obj=tanh(u),obj=sad_tree('tanh',{u}); end
        function obj=min(u,v),obj=sad_tree('min',{u,v}); end
        function obj=max(u,v),obj=sad_tree('max',{u,v}); end
        function obj=sum(u,v),obj=sad_tree('sum',{u,v}); end
        function obj=normpdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=sad_tree('normpdf',{u,v,w});
        end
        function obj=normcdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            obj=sad_tree('normcdf',{u,v,w});
        end
        function obj=abs(u),obj=sad_tree('abs',{u}); end
        function obj=isreal(u),obj=sad_tree('isreal',{u}); end
        % % % %         varagout=char(varargin)
        function string=char(obj,unravel,isparent)
            if nargin<3
                isparent=false;
                if nargin<2
                    unravel=false;
                end
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
                        operator=str2func(['my',obj.name]);
                        string=operator(mychar(args_{1},unravel),mychar(args_{2},unravel));
                    case 'uplus'
                        string=myplus('0',mychar(args_{1},unravel));
                    case 'uminus'
                        string=myminus('0',mychar(args_{1},unravel));
                    case {'mtimes'}
                        string=mytimes(mychar(args_{1},unravel),mychar(args_{2},unravel));
                    case {'mpower'}
                        string=mypower(mychar(args_{1},unravel),mychar(args_{2},unravel));
                    case {'rdivide','mrdivide'}
                        string=mydivide(mychar(args_{1},unravel),mychar(args_{2},unravel));
                    case {'min','max','gt','lt','ge','le'}
                        string=[obj.name,'(',mychar(args_{1},unravel),',',mychar(args_{2},unravel),')'];
                    case {'ldivide','mldivide'}
                        string=mydivide(mychar(args_{2},unravel),mychar(args_{1},unravel));
                    case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                            'tan','atan','tanh','abs','sqrt','isreal','sign'}
                        string=[obj.name,'(',mychar(args_{1},unravel),')'];
                    case {'normpdf','normcdf'}
                        string=[obj.name,'(',mychar(args_{1},unravel),',',mychar(args_{2},unravel),',',mychar(args_{3},unravel),')'];
                end
            end
            if isparent
                string=[obj.ref,'=',string];
            end
        end
        
        function obj=sqrt(u),obj=sad_tree('sqrt',{u}); end
        function obj=norm(u),obj=sad_tree('norm',{u}); end
        %         function obj=gt(u,v),obj=sad_tree('gt',{u,v}); end
        %         function obj=ge(u,v),obj=sad_tree('ge',{u,v}); end
        %         function obj=lt(u,v),obj=sad_tree('lt',{u,v}); end
        %         function obj=le(u,v),obj=sad_tree('le',{u,v}); end
        %         function obj=sign(u),obj=sad_tree('sign',{u}); end
        function this=diff(obj,wrt)
            if isnumeric(obj.name)
                this=sad_tree(0);
            elseif isempty(obj.args)
                if isequal(obj.name,wrt.name)
                    this=sad_tree(1);
                else
                    this=sad_tree(0);
                end
            else
                u=sad_tree(obj.args{1});
                du=diff(u,wrt);
                if numel(obj.args)>1
                    v=sad_tree(obj.args{2});
                    dv=diff(v,wrt);
                    if numel(obj.args)>2
                        w=sad_tree(obj.args{3});
                    end
                end
                switch obj.name
                    case 'plus'
                        this=du+dv;
                    case 'uplus'
                        this=du;
                    case 'minus'
                        this=du-dv;
                    case 'uminus'
                        this=-du;
                    case {'mtimes','times'}
                        update_reference(u);
                        update_reference(v);
                        upv=du*v;
                        vpu=dv*u;
                        this=upv+vpu;
                    case {'power','mpower'}
                        update_reference(u);
                        update_reference(v);
                        this=(log(u)+du/u*v)*u^v;
                    case {'rdivide','mrdivide'}
                        update_reference(u);
                        update_reference(v);
                        this=(du*v-dv*u)/v^2; % only when v~=0
                    case {'ldivide','mldivide'}
                        update_reference(u);
                        update_reference(v);
                        this=(u*dv-v*du)/u^2; % only when u~=0
                    case 'exp'
                        update_reference(obj);
                        this=du*obj;
                    case 'log'
                        update_reference(u);
                        this=du/u;
                    case 'log10'
                        update_reference(u);
                        this=(du/u)/log(10);
                    case 'cos'
                        update_reference(u);
                        this=-du*sin(u);
                    case 'acos'
                        update_reference(u);
                        this=-du/sqrt(1-u^2);
                    case 'cosh'
                        update_reference(u);
                        this=du*sinh(u);
                    case 'sin'
                        update_reference(u);
                        this=du*cos(u);
                    case 'asin'
                        update_reference(u);
                        this=du/sqrt(1-u^2);
                    case 'sinh'
                        update_reference(u);
                        this=du*cosh(u);
                    case 'tan'
                        update_reference(u);
                        this=du/(cos(u))^2;
                    case 'atan'
                        update_reference(u);
                        this=du/(1+u^2);
                    case 'tanh'
                        update_reference(u);
                        this=du/(1-u^2);
                    case 'min'
                        update_reference(u);
                        update_reference(v);
                        muv=u<v;
                        this=muv*u+(1-muv)*v;
                    case 'max'
                        update_reference(u);
                        update_reference(v);
                        muv=u>v;
                        this=muv*u+(1-muv)*v;
                    case 'sum'
                        this=sum(du);
                    case 'normpdf'
                        update_reference(u);
                        update_reference(v);
                        update_reference(w);
                        this=-du/w*(u-v)/w*obj;
                    case 'normcdf'
                        update_reference(u);
                        update_reference(v);
                        update_reference(w);
                        this=du*normpdf(u,v,w);
                    case 'abs' % this would not work without overloading sign!
                        update_reference(u);
                        this=du*sign(u);
                    case 'isreal'
                        update_reference(u);
                        this=isreal(u)*du;
                    case 'sqrt'
                        update_reference(obj);
                        this=du/(2*sqrt(obj));
                    case 'norm' % this would not work!
                        update_reference(u);
                        this=sum(u.*du)/norm(u);
                end
            end
        end
        function update_reference(obj)
            % not a variable and not a constant
            if isempty(obj.ref) && ~isempty(obj.args)
                obj.ncalls=obj.ncalls+1;
                if obj.ncalls>1
                    obj.ref='x';
                end
            end
        end
        function imax=re_flag_tree(tree,istart)
            if nargin<2
                istart=0;
            end
            if ~isempty(tree.args)
                for iarg=1:numel(tree.args)
                    istart=re_flag_tree(tree.args{iarg},istart);
                end
                if strcmp(tree.ref,'x')
                    istart=istart+1;
                    tree.ref=['T_',int2str(istart)];
                end
            end
            imax=istart;
        end
        function references=collect_references(tree,references)
            if nargin<2
                references=[];
            end
            if ~isempty(tree.args)
                isparent=~isempty(tree.ref) && is_atom(char(tree,false));
                if isparent
                    references=[{char(tree,false,true)};references];
                end
                for iarg=1:numel(tree.args)
                    references=collect_references(tree.args{iarg},references);
                end
            end
        end
        function print(tree)
            if ~isempty(tree.ref)
                fprintf(1,'%s\n',[tree.ref,' ---> ',tree.name]);
            end
            tree.args=reprocess_arguments(tree.args);
            for iarg=1:numel(tree.args)
                if ~ischar(tree.args{iarg})
                    print(tree.args{iarg})
                end
            end
        end
    end
    methods(Static)
        function varargout=re_sad_tree(varargin)
            varargout=varargin;
            for ii=1:length(varargin)
                if ~isa(vi,'sad_tree')
                    varargout{ii}=sad_tree(varargin{ii});
                end
            end
        end
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

function c=myplus(a,b)
if strcmp(a,'0')
    if strcmp(b,'0')
        c='0';
    else
        c=tryevaluate(b);
    end
else
    if strcmp(b,'0')
        c=tryevaluate(a);
    else
        c=tryevaluate([a,'+',b]);
    end
end
end

function c=myminus(a,b)
if strcmp(a,'0')
    if strcmp(b,'0')
        c='0';
    else
        c=tryevaluate(['-',parenthesize(b,'+-')]);
    end
else
    if strcmp(b,'0')
        c=a;
    else
        c=tryevaluate([a,'-',parenthesize(b,'+-')]);
    end
end
end

function c=mytimes(a,b)
if strcmp(a,'0')||strcmp(b,'0')
    c='0';
elseif strcmp(a,'1')
    c=b;
elseif strcmp(b,'1')
    c=a;
else
    c=tryevaluate([parenthesize(a,'+-'),'*',parenthesize(b,'+-')]);
end
end

function c=mypower(a,b)
if strcmp(b,'0')||strcmp(a,'1')
    c='1';
else
    c=tryevaluate([parenthesize(a,'+-/*^'),'^',parenthesize(b,'+-/*^')]);
end
end

function c=mydivide(a,b)
if strcmp(a,'0')
    c='0';
else
    c=tryevaluate([parenthesize(a,'+-'),'/',parenthesize(b,'+-/*^')]);
end
end

function x=parenthesize(x,forbid)
if nargin<2
    forbid='+-*/^';
end
forbid=strrep(forbid,'-','\-'); % must escape the minus sign
flag=~isempty(regexp(x,['[',forbid,']'],'start'));
if flag
    x=['(',x,')'];
end
end

function flag=is_atom(string)
flag=isempty(regexp(string,'[/*\-+^]','start'));
end

function a=tryevaluate(a)
% checks whether a string can be evaluated
flag=~any(isstrprop(a,'alpha'));
if flag
    cntrl=a;
    cntrl(isstrprop(cntrl,'digit'))=[];
    flag=~isempty(cntrl) && ~isequal(cntrl,'.');
    if flag
        flag=false;
        for ii=1:length(cntrl)
            if any(cntrl(ii)=='+-*^/')
                flag=true;
                break
            end
        end
        if flag
            a=num2str(eval(a),10);
        end
    end
end
end