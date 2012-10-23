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
            elseif ~unravel && ~isempty(obj.ref) && ~strcmp(obj.ref,'x') && ~isparent
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
                string=[obj.ref,'=',string,';'];
            end
        end
        function this=diff(obj,wrt)
            if isnumeric(obj.name)
                this=rise_sad(0);
            elseif isempty(obj.args)
                if isequal(obj.name,wrt.name)
                    this=rise_sad(1);
                else
                    this=rise_sad(0);
                end
            else
                u=rise_sad(obj.args{1});
                n_args=numel(obj.args);
                if n_args>1
                    v=rise_sad(obj.args{2});
                    if n_args>2
                        w=rise_sad(obj.args{3});
                    end
                end
                if iscell(wrt)
                    this=cell(size(wrt));
                    for ic=1:numel(wrt)
                        this{ic}=differentiation_engine(wrt{ic});
                    end
                elseif isa(wrt,'rise_sad')
                    this=differentiation_engine(wrt);
                else
                    error([mfilename,':: second argument must be rise_sad object'])
                end
            end
            function this=differentiation_engine(wrt)
                du=diff(u,wrt);
                switch obj.name
                    case {'gt','ge','lt','le','sign'}
                        this=rise_sad(0);
                    case 'plus'
                    dv=diff(v,wrt);
                        this=du+dv;
                    case 'uplus'
                        this=du;
                    case 'minus'
                    dv=diff(v,wrt);
                        this=du-dv;
                    case 'uminus'
                        this=-du;
                    case {'mtimes','times'}
                    dv=diff(v,wrt);
                        upv=du*v;
                        vpu=dv*u;
                        this=upv+vpu;
                        update_reference(u,v);
                    case {'power','mpower'}
                        this=(log(u)+du/u*v)*u^v;
                        update_reference(u,v);
                    case {'rdivide','mrdivide'}
                    dv=diff(v,wrt);
                        this=(du*v-dv*u)/v^2; % only when v~=0
                        update_reference(u,v);
                    case {'ldivide','mldivide'}
                    dv=diff(v,wrt);
                        this=(u*dv-v*du)/u^2; % only when u~=0
                        update_reference(u,v);
                    case 'exp'
                        this=du*obj;
                        update_reference(obj);
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
                        this=-du/w*(u-v)/w*obj;
                        update_reference(u,v,w);
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
                        this=du/(2*sqrt(obj));
                        update_reference(obj);
                    case 'norm' % this would not work!
                        this=sum(u.*du)/norm(u);
                        update_reference(u);
                    otherwise
                        keyboard
                end
            end
        end
        function update_reference(varargin)
            % not a variable and not a constant
            for iarg=1:length(varargin)
                if isempty(varargin{iarg}.ref) && ~isempty(varargin{iarg}.args)
                    varargin{iarg}.ncalls=varargin{iarg}.ncalls+1;
                    if varargin{iarg}.ncalls>1
                        varargin{iarg}.ref='x';
                    end
                end
            end
% % % % % % %             archive_it=~isnumeric(obj.name) && ~isempty(obj.args);
% % % % % % %             if archive_it
% % % % % % %                 if ~exist('mapObj_rise_sad','var')
% % % % % % %                     mapObj_rise_sad = containers.Map();
% % % % % % %                 end
% % % % % % %                 if ~isKey(mapObj_rise_sad,obj)
% % % % % % %                     nobj=size(mapObj_rise_sad,1);
% % % % % % %                     mapObj_rise_sad=containers.Map(obj,['T_',int2str(nobj)]);
% % % % % % %                 end
% % % % % % %                 assignin('base','mapObj_rise_sad',mapObj_rise_sad)
% % % % % % %             end
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
                isparent=~isempty(tree.ref) && strcmp(tree.ref(1),'T');% <---is_atom(char(tree,false));
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
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
        function varargout=re_sad_tree(varargin)
            varargout=varargin;
            for ii=1:length(varargin)
                if ~isa(vi,'rise_sad')
                    varargout{ii}=rise_sad(varargin{ii});
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