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
% % % %                 if ~exist('mapObj_rise_sad','var')
% % % %                     assignin('base','mapObj_rise_sad',containers.Map());
% % % %                     disp([mfilename,':: containers.map object created'])
% % % %                 end
% % % %                 mapObj_rise_sad=evalin('base','mapObj_rise_sad');
% % % %                 if ~isKey(mapObj_rise_sad,string)
% % % %                     mapObj_rise_sad=containers.Map(string,obj.ref);
% % % %                     assignin('base','mapObj_rise_sad',mapObj_rise_sad)
% % % %                 else
% % % %                     T_=values(mapObj_rise_sad,string);
% % % %                     if ~strcmp(T_,obj.ref)
% % % %                         obj.ref=T_;
% % % %                     end
% % % %                 end
                string=[obj.ref,'=',string,';']; % <---string={obj.ref,'=',string};
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
        function [this,references,tank,newtank]=diff(obj,wrt)
            % check that I do not need to output the object thanks to the
            % handle property
            [obj,tank]=consolidate(obj);
            this=mydiff(obj,wrt);
            newtank=containers.Map(values(tank),zeros(tank.Count,1));
            update_flags(obj,newtank);
            references=collect_references(obj);
        end
    end
    methods(Access=private)
        function update_flags(tree,newtank)
            % create new tank with all the T_ as keys and 0 as values for
            % all
            update_flags_intern(tree,1);
            update_flags_intern(tree,2);
            function update_flags_intern(tree,pass)
                for it=1:numel(tree)
                    if ~isempty(tree(it).args)
                        for iarg=1:numel(tree(it).args)
                            update_flags_intern(tree(it).args{iarg},pass);
                        end
                        calls=newtank(tree(it).ref);
                        if pass==1
                            if tree(it).ncalls
                                newtank(tree(it).ref)=calls+tree(it).ncalls;
                            end
                        else
                            if calls<2
                                tree(it).ref=[];
                            end
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
                        this=(log(u)+du/u*v)*u^v;
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
                        this=du/(2*sqrt(obj(index)));
                        update_reference(obj(index));
                    case 'norm' % this would not work!
                        this=sum(u.*du)/norm(u);
                        update_reference(u);
                    otherwise
                        keyboard
                end
            end
        end
        function [tree,tank]=consolidate(tree,tank)
            % vectorizes and consolidates a tree vector... 
            % I should check that I do not need to call the output because
            % of the handle property
            if nargin<2
                tank=containers.Map();
            end
            for it=1:numel(tree)
                for iarg=1:numel(tree(it).args)
                    [tree(it).args{iarg},tank]=consolidate(tree(it).args{iarg},tank);
                end
                if ~isempty(tree(it).args)
                    string=char(tree(it));
                    if isKey(tank,string)
                        T_=tank(string);
                    else
                        n=tank.Count+1;
                        T_=['T_',int2str(n)];
                        tank(string)=['T_',int2str(n)];
                    end
                    tree(it).ref=T_;
                end
            end
        end
        function update_reference(varargin)
            % not a variable and not a constant
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
            for it=1:numel(tree)
                if ~isempty(tree(it).args)
                    isparent=~isempty(tree(it).ref) && strcmp(tree(it).ref(1),'T');% <---is_atom(char(tree(it),false));
                    for iarg=1:numel(tree(it).args)
                        references=collect_references(tree(it).args{iarg},references);
                    end
                    if isparent
                        new_item=char(tree(it),false,true);
                        if ~any(strcmp(new_item,references))
                            references=[references;{new_item}]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
% % % %         function varargout=re_sad_tree(varargin)
% % % %             varargout=varargin;
% % % %             for ii=1:length(varargin)
% % % %                 if ~isa(vi,'rise_sad')
% % % %                     varargout{ii}=rise_sad(varargin{ii});
% % % %                 end
% % % %             end
% % % %         end
% % % %         function [tree,tank]=match_tree(tree,tank)
% % % %             for it=1:numel(tree)
% % % %                 nail=tree(it)==tank;
% % % %                 if any(nail)
% % % %                     tree(it)=tank(nail);
% % % %                     disp('match found and successfully pushed')
% % % %                 else
% % % %                     % grow the tank
% % % %                     tank=[tank;tree(it)]; %#ok<AGROW>
% % % %                     nargs=numel(tree(it).args);
% % % %                     for iarg=1:nargs
% % % %                         [tree(it).args{iarg},tank]=match_tree(tree(it).args{iarg},tank);
% % % %                     end
% % % %                 end
% % % %             end
% % % %         end
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

%{
clear all
clc
a=rise_sad('a');b=rise_sad('b');c=rise_sad('c');
func=@(a,b,c)[exp(a+2*log(b+c)-a*atan(b*c));exp(-a*atan(b*c));exp(a+2*log(b+c))];
tree=func(a,b,c);
tmp=tree.vectorize
tmp==tmp(2)
*while vectorizing, check equalities?
* construct a new tree and check that its elements are not in the previous
trees for consistency and savings... then 
dd=diff(tree,{a,b,c});
re_flag_tree(tree)
references=collect_references(tree)
%}
