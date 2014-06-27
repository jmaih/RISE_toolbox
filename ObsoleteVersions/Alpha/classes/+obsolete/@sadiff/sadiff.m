classdef sadiff % < handle
    % to do: unknown functions
    properties
        func
        % name of variable or function
        args
        % cell array containing the arguments
        order
        % order of variable in the list of the vars to be differentiated
        % wrt. This is useful especially when computing higher-order
        % derivatives, in order to avoid computing the same derivative
        % twice.
    end
    methods
        function obj=sadiff(func_,args_,order_)
            if nargin
                if isa(func_,'sadiff')
                    obj=func_;
                    return
                end
                if nargin<3
                    order_=[];
                    if nargin<2
                        args_=cell(0);
                    end
                    if ~iscell(args_)
                        args_={args_};
                    end
                end
                obj.func=func_;
                obj.args=args_;
                obj.func=func_;
                obj.order=order_;
            end
        end
        function obj=plus(u,v)
            [a,b]=sadiff.commute(u,v);
            obj=sadiff.prototype('+',a,b); 
        end
        function obj=uplus(u),obj=sadiff.prototype('uplus',u); end
        function obj=minus(u,v),obj=sadiff.prototype('-',u,v); end
        function obj=uminus(u),obj=sadiff.prototype('uminus',u); end
        function obj=times(u,v),obj=mtimes(u,v); end
        function obj=mtimes(u,v)
            [a,b]=sadiff.commute(u,v);
            obj=sadiff.prototype('*',a,b); 
        end
        function obj=power(u,v),obj=mpower(u,v); end
        function obj=mpower(u,v),obj=sadiff.prototype('^',u,v); end
        function obj=rdivide(u,v),obj=mrdivide(u,v); end
        function obj=mrdivide(u,v),obj=sadiff.prototype('/',u,v); end
        function obj=ldivide(u,v),obj=mldivide(u,v); end
        function obj=mldivide(u,v),obj=mrdivide(v,u); end
        function obj=exp(u),obj=sadiff.prototype('exp',u); end
        function obj=log(u),obj=sadiff.prototype('log',u); end
        function obj=log10(u),obj=log(u)/log(10); end
        function obj=cos(u),obj=sadiff.prototype('cos',u); end
        function obj=acos(u),obj=sadiff.prototype('acos',u); end
        function obj=cosh(u),obj=sadiff.prototype('cosh',u); end
        function obj=sin(u),obj=sadiff.prototype('sin',u); end
        function obj=asin(u),obj=sadiff.prototype('asin',u); end
        function obj=sinh(u),obj=sadiff.prototype('sinh',u); end
        function obj=tan(u),obj=sadiff.prototype('tan',u); end
        function obj=atan(u),obj=sadiff.prototype('atan',u); end
        function obj=tanh(u),obj=sadiff.prototype('tanh',u); end
        function obj=sqrt(u),obj=sadiff.prototype('sqrt',u); end
        function obj=normpdf(u,mu,sig)
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            obj=sadiff.prototype('normpdf',u,mu,sig);
        end
        function obj=normcdf(u,mu,sig)
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            obj=sadiff.prototype('normcdf',u,mu,sig);
        end
        function names=load_varlist(obj,names)
            if nargin<2
                names={};
            end
            if isa(obj,'sadiff')
                if isempty(obj.args)
                    if ischar(obj.func)
                        names=union(names,obj.func);
                    end
                else
                    for iarg=1:numel(obj.args)
                        if isa(obj.args{iarg},'sadiff')
                            names=load_varlist(obj.args{iarg},names);
                        end
                    end
                end
            end
        end
         function this=expand(obj)
            this=sadiff.empty(0);
            for iobj=1:size(obj,1)
                this(iobj,:)=re_expand(obj(iobj));
            end
             function this=re_expand(this)
                this_args=this.args;
                nargs=numel(this_args);
                if nargs
                    max_args=0;
                    for iarg=1:nargs
                        switch class(this_args{iarg})
                            case 'sadiff'
                                this_args{iarg}=re_expand(this_args{iarg});
                                this_args{iarg}=num2cell(this_args{iarg});
                                max_args=max(max_args,numel(this_args{iarg}));
                            case {'double','char'}
                                this_args{iarg}=this_args(iarg);
                        end
                    end
                    max_args=max(max_args,1);
                    s=cell(max_args,nargs);
                    for iarg=1:nargs
                        arg_this=this_args{iarg};
                        if numel(arg_this)<max_args
                            arg_this=arg_this(ones(1,max_args));
                        end
                        s(:,iarg)=arg_this(:);
                    end
                    this=this(ones(1,max_args));
                    % this is where we should remove 0 and 1 is applicable.
                    % that is 0* 0/ +0 -0 ^0 *1 /1 ^1
                    for ii=1:max_args
                        this(ii).args=s(ii,:);
                    end
                end
            end
        end
        varargout=differentiate(varargin)
        varargout=diff(varargin)
        varargout=char(varargin)
        varargout=print(varargin)
   end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
        varargout=setup(varargin)
    end
    methods(Static,Access=private)
        varargout=neat(varargin)
        function [a,b]=commute(x,y)
            a=x;
            b=y;
            if ~isa(a,'sadiff') && ~isa(b,'sadiff')
                testa=get_tester(a);
                testb=get_tester(b);
                swap=testa>testb;
                if swap
                    a=y;
                    b=x;
                end
            end
            function zz=get_tester(z)
                if ischar(z)
                    zz=z(1);
                elseif isa(z,'double')
                    zz=z;
                else
                    error('unknown type')
                end
            end
        end
        function obj=prototype(func,varargin)
            persistent modell__
            if isempty(modell__)
                modell__=sadiff('xxx');
            end
            obj=modell__;
            if nargin
                obj.func=func;
                if nargin>1
                    obj.args=varargin;
                end
            end
        end
        function mycall=metastruct()
            mycall=struct('operCount',0,'lineCount',0,'fid',{cell(0,2)},'prefix_list',{{'xx','dF','F','indx'}});
        end
        function mycall=trim_metastruct(mycall,restricted_list)
            oldcall=mycall;
            while 1
                mycall=trim_engine(mycall,restricted_list);
                if isequal(oldcall,mycall)
                    break
                end
                oldcall=mycall;
            end
            function mycall=trim_engine(mycall,restricted_list)
                prefix=cell2mat(strcat(mycall.prefix_list,'|'));
                test=regexp(mycall.fid(:,2),['(?<!\w)(',prefix(1:end-1),')_[0-9]+_(?!\w)'],'match');
                test=[test{:}];
                lineCount=size(mycall.fid,1);
                discard=false(lineCount,1);
                for iop=1:lineCount
                    xi=mycall.fid{iop,1};
                    if any(strcmp(xi,restricted_list))
                        continue
                    end
                    weight=sum(strcmp(xi,test));
                    if weight<1 % discard if the variable does not appear anywhere on the rhs
                        discard(iop)=true;
                    end
                end
                mycall.fid=mycall.fid(~discard,:);
                mycall.lineCount=size(mycall.fid,1);
            end
        end
    end
end

