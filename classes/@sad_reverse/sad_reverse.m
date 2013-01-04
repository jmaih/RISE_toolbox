classdef sad_reverse
    properties
        name
        args
        order
    end
    methods
        %------------------
        % constructor
        %------------------
        function obj=sad_reverse(name,arg,order)
            if nargin
                if isa(name,'sad_reverse')
                    obj=name;
                else
                    if nargin<3
                        order=[];
                        if nargin<2
                            arg=[];
                        end
                    end
                    if isempty(arg),arg={}; end
                    if ~iscell(arg),arg={arg}; end
                    obj.name=name;
                    obj.args=arg;
                    obj.order=order;
                end
            end
        end
        %------------------
        % unary functions
        %------------------
        function obj=abs(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@abs,u);
                return
            end
            obj=sad_reverse('abs',{u});
        end
        function obj=isreal(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@isreal,u);
                return
            end
            obj=sad_reverse('isreal',{u});
        end
        function obj=sqrt(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@sqrt,u);
                return
            end
            if isnumeric(u)
                obj=sad_reverse(sqrt(u.name));
            else
                obj=sad_reverse('sqrt',{u});
            end
        end
        function obj=sign(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@sign,u);
                return
            end
            obj=sad_reverse('sign',{u});
        end
        function obj=uplus(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@uplus,u);
                return
            end
            obj=u;
        end
        function obj=uminus(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@uminus,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(0);
                %test1=dbstack;disp([%test1(1).name,' :: zu'])
            elseif isequal(u.name,'uminus')
                obj=sad_reverse(u.args{1}); % -(-u)=u
                %test1=dbstack;disp([%test1(1).name,' :: -(-u)'])
            else
                obj=sad_reverse('uminus',{u});
            end
        end
        function obj=exp(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@exp,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(1);
            elseif isone(u)
                obj=sad_reverse(exp(1));
            else
                obj=sad_reverse('exp',{u});
            end
        end
        function obj=log(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@log,u);
                return
            end
            if isone(u)
                obj=sad_reverse(0);
            elseif iszero(u)
                error('log of zero')
            else
                obj=sad_reverse('log',{u});
            end
        end
        function obj=log10(u)
            obj=u/log(10);
        end
        function obj=cos(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@cos,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(1);
            else
                obj=sad_reverse('cos',{u});
            end
        end
        function obj=acos(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@acos,u);
                return
            end
            if isone(u)
                obj=sad_reverse(0);
            else
                obj=sad_reverse('acos',{u});
            end
        end
        function obj=cosh(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@cosh,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(1);
            else
                obj=sad_reverse('cosh',{u});
            end
        end
        function obj=sin(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@sin,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(0);
            else
                obj=sad_reverse('sin',{u});
            end
        end
        function obj=asin(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@asin,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(1);
            else
                obj=sad_reverse('asin',{u});
            end
        end
        function obj=sinh(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@sinh,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(0);
            else
                obj=sad_reverse('sinh',{u});
            end
        end
        function obj=tan(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@tan,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(0);
            else
                obj=sad_reverse('tan',{u});
            end
        end
        function obj=atan(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@atan,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(0);
            else
                obj=sad_reverse('atan',{u});
            end
        end
        function obj=tanh(u)
            u=sad_reverse(u);
            if numel(u)>1
                obj=sad_reverse.recursive(@tanh,u);
                return
            end
            if iszero(u)
                obj=sad_reverse(0);
            else
                obj=sad_reverse('tanh',{u});
            end
        end
        %------------------
        % binary functions
        %------------------
        function obj=plus(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@plus,u,v);
                return
            end
            % try a quick exit
            if ~isnumeric(v) && strcmp(v.name,'uminus')
                obj=minus(u,sad_reverse(v.args{1}));
            else
                zu=iszero(u); zv=iszero(v);
                if zu && zv
                    obj=u;
                elseif isequal(u,v)
                    obj=2*u;
                elseif zu
                    obj=v;
                elseif zv
                    obj=u;
                else
                    obj=sad_reverse('plus',{u,v});
                end
            end
        end
        function obj=minus(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@minus,u,v);
                return
            end
            % try a quick exit
            zu=iszero(u); zv=iszero(v);
            if (zu && zv)||isequal(u,v)
                obj=sad_reverse(0);
            elseif zu
                obj=-v;
            elseif zv
                obj=u;
            else
                obj=sad_reverse('minus',{u,v});
            end
        end
        function obj=times(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@times,u,v);
                return
            end
            % try a quick exit
            zu=iszero(u); zv=iszero(v);
            if zu || zv
                obj=sad_reverse(0);
            elseif isequal(u,v)
                obj=u^2;
            else
                if isone(v)
                    obj=u;
                elseif isone(u)
                    obj=v;
                else
                    obj=sad_reverse('times',{u,v});
                end
            end
        end
        function obj=mtimes(u,v)
            obj=times(u,v);
        end
        function obj=power(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@power,u,v);
                return
            end
            % try a quick exit
            if iszero(v) || isone(u)
                obj=sad_reverse(1);
                %test1=dbstack;disp([%test1(1).name,' :: v=0||u=1'])
            elseif isone(v)
                obj=sad_reverse(u);
                %test1=dbstack;disp([%test1(1).name,' :: v=1'])
            else
                obj=sad_reverse('power',{u,v});
            end
        end
        function obj=mpower(u,v)
            obj=power(u,v);
        end
        function obj=rdivide(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@rdivide,u,v);
                return
            end
            % try a quick exit
            if isone(v)
                obj=u;
                %test1=dbstack;disp([%test1(1).name,' :: v=1'])
            elseif iszero(v)
                error('division by zero')
            elseif iszero(u)
                obj=u;
                %test1=dbstack;disp([%test1(1).name,' :: u=0'])
            elseif isequal(u,v)
                obj=sad_reverse(1);
                %test1=dbstack;disp([%test1(1).name,' :: u=v'])
            else
                obj=sad_reverse('rdivide',{u,v});
            end
        end
        function obj=mrdivide(u,v)
            obj=rdivide(u,v);
        end
        function obj=ldivide(u,v)
            obj=rdivide(v,u);
        end
        function obj=mldivide(u,v)
            obj=ldivide(u,v);
        end
        function obj=min(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@min,u,v);
                return
            end
            obj=sad_reverse('min',{u,v});
        end
        function obj=max(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@max,u,v);
                return
            end
            obj=sad_reverse('max',{u,v});
        end
        function obj=gt(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@gt,u,v);
                return
            end
            obj=sad_reverse('gt',{u,v});
        end
        function obj=ge(u,v)
            % try a quick exit
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@ge,u,v);
                return
            end
            obj=sad_reverse('ge',{u,v});
        end
        function obj=lt(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@lt,u,v);
                return
            end
            obj=sad_reverse('lt',{u,v});
        end
        function obj=le(u,v)
            v=sad_reverse(v);
            u=sad_reverse(u);
            if numel(u)>1||numel(v)>1
                obj=sad_reverse.recursive(@le,u,v);
                return
            end
            obj=sad_reverse('le',{u,v});
        end
        %------------------
        % trinary functions
        %------------------
        function obj=normpdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            u=sad_reverse(u);
            v=sad_reverse(v);
            w=sad_reverse(w);
            if numel(u)>1||numel(v)>1||numel(w)>1
                obj=sad_reverse.recursive(@normpdf,u,v,w);
                return
            end
            obj=sad_reverse('normpdf',{u,v,w});
        end
        function obj=normcdf(u,v,w)
            if nargin<3
                w=1;
                if nargin<2
                    v=0;
                end
            end
            u=sad_reverse(u);
            v=sad_reverse(v);
            w=sad_reverse(w);
            if numel(u)>1||numel(v)>1||numel(w)>1
                obj=sad_reverse.recursive(@normcdf,u,v,w);
                return
            end
            obj=sad_reverse('normcdf',{u,v},w);
        end
        %------------------
        % helper functions
        %------------------
        function flag=isnumeric(u)
            flag=nan(size(u));
            for ii=1:numel(u)
                flag(ii)=isnumeric(u(ii).name);
            end
        end
        function flag=iszero(u)
            sizu=size(u);
            flag=zeros(sizu);
            for ii=1:prod(sizu)
                flag(ii)=isequal(u(ii).name,0);
            end
        end
        function flag=isone(u)
            flag=nan(size(u));
            for ii=1:numel(u)
                flag(ii)=isequal(u(ii).name,1);
            end
        end
        %------------------
        % printing functions
        %------------------
        function operCount=print(tree)
            if nargin<2
                operCount=0;
            end
            for it=1:numel(tree)
                fprintf(1,'%s %0.16g\n','object number ',it);
                operCount=operCount+reprint(tree(it));
            end
            
            function operCount=reprint(obj,indent,last,operCount)
                if nargin<4
                    operCount=0;
                    if nargin<3
                        last=true;
                        if nargin<2
                            indent='';
                        end
                    end
                end
                if isnumeric(obj)
                    % do nothing?
                    obj.name=sprintf('%0.16g',obj.name);
                end
                nn='';
                if ~isempty(obj.args)
                    operCount=operCount+1;
                    nn='node::';
                end
                memo=indent;
                if last
                    memo=[memo,'\'];
                    indent =[indent,' '];
                else
                    memo=[memo,'|'];
                    indent =[indent,'| '];
                end
                disp([memo,nn,obj.name]);
                children_count=numel(obj.args);
                indentChild=indent; %<--- indentChild=[' ',indent];
                for iarg=1:children_count
                    obj.args{iarg}=sad_reverse(obj.args{iarg});
                    operCount=reprint(obj.args{iarg},indentChild,iarg==children_count,operCount);
                end
            end
        end
        function [c,map]=char(u,map,first_call)
            if nargin<3
                first_call=true;
                if nargin<2
                    map=[];
                end
            end
            if first_call
                c=cell(size(u));
                first_call=false;
                nobj=numel(u);
                for iobj=1:nobj
                    if ~isempty(map) && ~isempty(map.index) && iobj<=size(map.index,1)
                        % if the form is  not vectorized, the index will be
                        % short and will eventually be trimmed away anyway.
                        % Ideally it should not be included in the first
                        % place.
                        right=mat2str(map.index{iobj,2});
                        right(isspace(right))=',';
                        map=sad_reverse.update_map(map,right,map.index{iobj,1});
                    end
                    [c{iobj},map]=char(u(iobj),map,first_call);
                end
                if ~isempty(map)
                    map=sad_reverse.trim_metastruct(map,c);
                end
            else
                is_map_update=true;
                if isnumeric(u)
                    c=sprintf('%0.16g',u.name);
                    is_map_update=false;
                elseif isempty(u.args)
                    c=u.name;
                    % this line should be removed but will stay around for
                    % a while at least until I find a workaround for the
                    % case where there are not many derivatives to
                    % compute...
                    is_map_update=false; 
                else
                    u_args=u.args;
                    for iarg=1:numel(u_args)
                        switch class(u_args{iarg})
                            case 'char'
                                % do nothing that is where we are going and you
                                % are already there
                            case 'sad_reverse'
                                [u_args{iarg},map]=char(u_args{iarg},map,first_call);
                            case 'double'
                                u_args{iarg}=sprintf('%0.16g',u_args{iarg});
                            otherwise
                                error([class(u_args{iarg}),' unexpected at this stage. please contact junior.maih@gmail.com'])
                        end
                    end
                    c=sad_reverse.neat(u.name,u_args{:});
                    % a formal operation has taken place. We can update the
                    % map.
                end
                if is_map_update
                    [map,c]=sad_reverse.update_map(map,c);
                end
            end
        end
        function list=get_varlist(obj,list)
            if nargin<2
                list={};
            end
            if numel(obj)>1
                error('only one object at a time')
            end
            if isempty(obj.args)
                if ischar(obj.name)
                    list=union(list,obj.name);
                end
            else
                for iarg=1:numel(obj.args)
                    list=get_varlist(obj.args{iarg},list);
                end
            end
        end
        %------------------
        % differentiator
        %------------------
        function [this,index]=diff(obj,wrt,vectorize)
            if nargin<3
                vectorize=false;
            end
            % check that I do not need to output the object thanks to the
            % handle property
            if ischar(wrt)
                wrt=cellstr(wrt);
            end
            if iscellstr(wrt)
                wrt=sad_reverse.symbols(wrt,true);
            end
            if iscell(wrt)
                wrt=[wrt{:}];
            end
            if ~isa(wrt,'sad_reverse')
                error([mfilename,':: second argument must be sad_reverse object, a cellstr or a char'])
            end
            nwrt=numel(wrt);
            for irt=1:nwrt
                if vectorize && isempty(wrt(irt).order)
                    error('the order of the variables should be provided for vectorized differentiation')
                end
            end
            orders=1:nwrt;
            if vectorize && ~isequal(orders,[wrt.order])
                error('something wrong: somehow list of wrt does not match wrt')
            end
            wrt=wrt(:)';
            wrt_names={wrt.name};
            obj=obj(:);
            nobj=numel(obj);
            if vectorize
                index=cell(nobj,2);
                this=sad_reverse.empty(0,1);
                tmp=sad_reverse.create_map();
                prefix=tmp.prefix_list{2};
                clear tmp
            else
                index=[];
                this=sad_reverse(0);
                this=this(ones(nobj,nwrt));
            end
            for iobj=1:nobj
                mylist=get_varlist(obj(iobj));
                tt=ismember(wrt_names,mylist);
                targets=orders(tt);
                if isempty(targets)
                    continue
                end
                if vectorize
                    index(iobj,:)={[prefix,'_',int2str(iobj),'_'],...
                        targets};
                    this(iobj,1)=mydiff(obj(iobj),wrt(targets),index{iobj,1});
                else
                    this(iobj,targets)=mydiff(obj(iobj),wrt(targets));
                end
            end
        end
    end
    %----------------------
    % Private methods
    %----------------------
    methods(Access=private)
        function this=mydiff(obj,wrt,indx_name)
            if nargin<3
                indx_name=[];
            end
            vectorize=~isempty(indx_name);
            if isnumeric(obj.name)
                this=sad_reverse(0);
                % derivatives =0 as above
            elseif isempty(obj.args)
                nwrt=numel(wrt);
                this=sad_reverse(0);
                if vectorize
                    loc=find(strcmp(obj.name,{wrt.name}));
                    if ~isempty(loc)
                        orders=[wrt.order];
                        left=sprintf('%0.16g',orders(loc));
                        this=sad_reverse(['bigi_(',left,',',indx_name,')']);
                    end
                else
                    this=this(ones([1,nwrt]));
                    wrt_names={wrt.name};
                    this(strcmp(obj.name,wrt_names))=sad_reverse(1);
                end
            else
                this=differentiation_engine(obj,wrt,indx_name);
            end
        end
    end
    %----------------------
    % Static private methods
    %----------------------
    methods(Static,Access=private)
        function this=recursive(func,varargin)
            nvars=length(varargin);
            siz_args=ones(1,nvars);
            cellarray=cell(1,nvars);
            state=ones(1,nvars);
            for iarg=1:nvars
                siz_args(iarg)=numel(varargin{iarg});
                cellarray{iarg}=varargin{iarg}(1);
            end
            this=func(cellarray{:});
            %             this=sad_reverse(func,cellarray);
            %             this=repmat(this,1,nvars);
            for iter=2:max(siz_args)
                state_iter=state(1,:);
                for iarg=1:nvars
                    if siz_args(iarg)>1
                        cellarray{iarg}=varargin{iarg}(iter);
                        state_iter(iarg)=iter;
                    end
                end
                test=sum(bsxfun(@minus,state,state_iter),2);
                loc=find(test==0);
                if isempty(loc)
                    this(1,iter)=func(cellarray{:});
                    state=[state;state_iter];
                else
                    this(1,iter)=this(1,loc);
                end
                %this(1,iter)=this(1).set_options('args',cellarray);
                %                 this(1,iter)=sad_reverse(func,cellaray);
                %                 this(iter).args=cellaray;
            end
        end
        function mycall=trim_metastruct(mycall,c,cutoff)
            if nargin<3
                cutoff=2;
            end
            restricted_list=sad_reverse.extract_list(c,mycall.prefix_list);
            % add the indx_ to the restricted list
            restricted_list=union(restricted_list,...
                sad_reverse.extract_list(mycall.fid(:,1),mycall.prefix_list(end)));
            % add some potential previously fixed references
            restricted_list=union(restricted_list,mycall.fixed_references);
            
            oldcall=mycall;
            while 1
                mycall=trim_engine(mycall,restricted_list);
                if isequal(oldcall,mycall)
                    break
                end
                oldcall=mycall;
            end
            function mycall=trim_engine(mycall,restricted_list)
                discard_list=sad_reverse.extract_list(mycall.fid(:,2),mycall.prefix_list);
                line_count=mycall.line_count;
                weights=zeros(line_count,1);
                for iop=1:line_count
                    xi=mycall.fid{iop,1};
                    if any(strcmp(xi,restricted_list))
                        weights(iop)=inf;
                    else
                        weights(iop)=sum(strcmp(xi,discard_list));
                    end
                end
                for iop=line_count:-1:1
                    if weights(iop)<=cutoff
                        oldstr=mycall.fid{iop,1};
                        newstr=mycall.fid{iop,2};
                        if sad_reverse.check_offend(newstr,'+-*/^');
                            newstr=['(',newstr,')']; %#ok<AGROW>
                        end
                        mycall.fid=strrep(mycall.fid,oldstr,newstr);
                        % some references like indx_ appear both on the
                        % left and on the right and so I have to replace
                        % both sides instead of just the rhs
                    end
                end
                discard=weights<=cutoff;
                mycall.fid=mycall.fid(~discard,:);
                mycall.line_count=size(mycall.fid,1);
            end
        end
        function list=extract_list(object,names)
            prefix=cell2mat(strcat(names,'|'));
            list=regexp(object,...
                ['(?<!\w)(',prefix(1:end-1),')_[0-9]+_(?!\w)'],'match');
            list=[list{:}];
        end
        function c=neat(func,varargin)
            fragile=true;
            switch func
                case {'abs','cos','sqrt','sin','acos','cosh',...
                        'tan','tanh','log10','isreal','asin','atan','log',...
                        'exp','uplus','sign','sinh'}
                    c=[func,'(',varargin{1},')'];
                case 'uminus'
                    if fragile==true && sad_reverse.check_offend(varargin{1},'+-')
                        varargin{1}=['(',varargin{1},')'];
                    end
                    c=['-',varargin{1}];
                case {'+','plus'}
                    if fragile==true
                        if (strcmp(varargin{1},'0') && strcmp(varargin{2},'0'))
                            c='0';
                        elseif strcmp(varargin{1},'0')
                            c=varargin{2};
                        elseif strcmp(varargin{1},varargin{2})
                            if sad_reverse.check_offend(varargin{1},'+-')
                                varargin{1}=['(',varargin{1},')'];
                            end
                            c=['2*',varargin{1}];
                        elseif strcmp(varargin{2},'0')
                            c=varargin{1};
                        else
                            c=[varargin{1},'+',varargin{2}];
                        end
                    else
                        c=[varargin{1},'+',varargin{2}];
                    end
                case {'-','minus'}
                    if fragile==true
                        if strcmp(varargin{1},varargin{2})
                            c='0';
                        elseif strcmp(varargin{1},'0')
                            if sad_reverse.check_offend(varargin{2},'+-')
                                varargin{2}=['(',varargin{2},')'];
                            end
                            c=['-',varargin{2}];
                        elseif strcmp(varargin{2},'0')
                            c=varargin{1};
                        else
                            if sad_reverse.check_offend(varargin{2},'+-')
                                varargin{2}=['(',varargin{2},')'];
                            end
                            c=[varargin{1},'-',varargin{2}];
                        end
                    else
                        c=[varargin{1},'-',varargin{2}];
                    end
                case {'*','times','mtimes'}
                    if fragile==true
                        if strcmp(varargin{1},'0') || strcmp(varargin{2},'0')
                            c='0';
                        elseif strcmp(varargin{1},'1') && strcmp(varargin{2},'1')
                            c='1';
                        elseif strcmp(varargin{1},'1')
                            c=varargin{2};
                        elseif strcmp(varargin{2},'1')
                            c=varargin{1};
                        elseif strcmp(varargin{1},varargin{2})
                            if really_neat  && sad_reverse.check_offend(varargin{1},'+-/*^')
                                varargin{1}=['(',varargin{1},')'];
                            end
                            c=[varargin{1},'.^2'];
                        else
                            if sad_reverse.check_offend(varargin{1},'+-')
                                varargin{1}=['(',varargin{1},')'];
                            end
                            if sad_reverse.check_offend(varargin{2},'+-')
                                varargin{2}=['(',varargin{2},')'];
                            end
                            c=[varargin{1},'.*',varargin{2}];
                        end
                    else
                        c=[varargin{1},'.*',varargin{2}];
                    end
                case {'/','rdivide','mrdivide'}
                    if fragile==true
                        if strcmp(varargin{1},'0')
                            c='0';
                        elseif strcmp(varargin{1},varargin{2})
                            c='1';
                        elseif strcmp(varargin{2},'1')
                            c=varargin{1};
                        else
                            if sad_reverse.check_offend(varargin{1},'+-')
                                varargin{1}=['(',varargin{1},')'];
                            end
                            if sad_reverse.check_offend(varargin{2},'+-*/')
                                varargin{2}=['(',varargin{2},')'];
                            end
                            c=[varargin{1},'./',varargin{2}];
                        end
                    else
                        c=[varargin{1},'./',varargin{2}];
                    end
                case {'^','power','mpower'}
                    if fragile==true
                        if strcmp(varargin{2},'0')
                            c='1';
                        elseif strcmp(varargin{2},'1')||strcmp(varargin{1},'1')
                            c=varargin{1};
                        else
                            if sad_reverse.check_offend(varargin{1},'+-*/^')
                                varargin{1}=['(',varargin{1},')'];
                            end
                            if sad_reverse.check_offend(varargin{2},'+-*/^')
                                varargin{2}=['(',varargin{2},')'];
                            end
                            c=[varargin{1},'.^',varargin{2}];
                        end
                    else
                        c=[varargin{1},'.^',varargin{2}];
                    end
                case {'normpdf','normcdf'}
                    c=[func,'(',varargin{1},',',varargin{2},',',varargin{3},')'];
                case {'\','ldivide','mldivide'}
                    if fragile==true
                        c=['(',varargin{2},')./(',varargin{1},')'];
                    else
                        c=[varargin{2},'./',varargin{1}];
                    end
                case {'lt','<'}
                    c=[varargin{1},'<',varargin{2}];
                case {'gt','>'}
                    c=[varargin{1},'>',varargin{2}];
                case {'le','<='}
                    c=[varargin{1},'<=',varargin{2}];
                case {'ge','>='}
                    c=[varargin{1},'>=',varargin{2}];
                case {'min','max'}
                    c=[func,'(',varargin{1},',',varargin{2},')'];
                otherwise
                    c=cell2mat(strcat(varargin,','));
                    c=[func,'(',c(1:end-1),')'];
            end
        end
        function [map,c]=update_map(map,c,imposed_name)
            if nargin<3
                imposed_name=[];
            end
            if isempty(map)
                return
            end
            prefix=map.prefix_list{1};
            if iscell(c)
                for icell=1:numel(c)
                    [map,c{icell}]=update_map(map,c{icell});
                end
            else
                loc=find(strcmp(c,map.fid(:,2)));
                if isempty(loc)
                    map.operation_count=map.operation_count+1;
                    map.line_count=map.line_count+1;
                    if ~isempty(imposed_name)
                        ref=imposed_name;
                    else
                        ref=[prefix,'_',sprintf('%0.16g',map.operation_count),'_'];
                    end
                    map.fid=[map.fid;{ref,c}];
                    c=ref;
                else
                    c=map.fid{loc,1};
                end
            end
        end
    end
    %----------------------
    % Static public methods
    %----------------------
    methods(Static)
        function flag=check_offend(v,offend)
            flag=false;
            iter=0;
            depth=0;
            while iter<length(v)
                iter=iter+1;
                if strcmp(v(iter),'(')
                    depth=depth+1;
                elseif strcmp(v(iter),')')
                    depth=depth-1;
                end
                if depth==0 && any(v(iter)==offend)
                    flag=true;
                    break
                end
            end
            % flag=any(ismember(v,offend));
        end
        function cellarray=symbols(cellarray,ordered)
            if nargin<2
                ordered=false;
            end
            if ischar(cellarray)
                cellarray=cellstr(cellarray);
            end
            if ~iscellstr(cellarray)
                error('elements must be char or cellstr');
            end
            for ic=1:numel(cellarray)
                cellarray{ic}=sad_reverse(cellarray{ic});
                if ordered
                    cellarray{ic}.order=ic;
                end
            end
        end
        function map=create_map(index)
            map=struct('fid',{cell(0,2)},'operation_count',0,'index',[],...
                'prefix_list',{{'ref','indx'}},'line_count',0,...
                'fixed_references',{{}}); % reference holders for higher-order derivatives
            if nargin
                map.index=index;
            end
        end
        function mat=create_matrix(map,c,nwrt,name)
            if ischar(c)
                c=cellstr(c);
            end
            [nrows,ncols]=size(c);
            prologue={name,['zeros(',sprintf('%0.16g',nrows),...
                ',',sprintf('%0.16g',nwrt),')']};
            is_vectorized=~isempty(map.index);
            if is_vectorized
                if ncols~=1
                    error('more than one column is unexpected. Contact junior.maih@gmail.com')
                end
                prologue=[prologue;{'bigi_',['speye(',sprintf('%0.16g',nwrt),')']}];
            else
                if nwrt<ncols
                    error('number of wrt variables < # cols of argument 2')
                end
                cols_str=cell(1,ncols);
            end
            epilogue=map.fid;
            for irow=1:nrows
                row_str=sprintf('%0.16g',irow);
                if is_vectorized
                    right=[map.prefix_list{2},'_',row_str,'_'];
                    str=[name,'(',row_str,',',right,')'];
                    % replace only the lhs of the epilogue
                    epilogue(:,1)=strrep(epilogue(:,1),c{irow},str);
                else
                    for icol=1:ncols
                        if irow==1
                            cols_str{icol}=sprintf('%0.16g',icol);
                        end
                        if ~isequal(c{irow,icol},'0')
                            str=[name,'(',row_str,',',cols_str{icol},')'];
                            %                             if strncmp('ref_',c{irow,icol},4)
                            %                                 % replace only the lhs of the epilogue
                            %                             epilogue(:,1)=strrep(epilogue(:,1),c{irow,icol},str);
                            %                             else
                            epilogue=[epilogue;{str,c{irow,icol}}]; %#ok<AGROW>
                            %                             end
                        end
                    end
                end
            end
            epilogue=[prologue;epilogue];
            mat=strcat(epilogue(:,1),'=',epilogue(:,2),';');
        end
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
    end
end

function this=differentiation_engine(myobj,wrt,indx_name)
args_=myobj.args;
for iarg=1:numel(args_)
    args_{iarg}=sad_reverse(args_{iarg});
end
u=args_{1};
du=mydiff(u,wrt,indx_name);
switch myobj.name
    case {'gt','ge','lt','le','sign'}
        this=sad_reverse(0);
    case 'plus'
        v=args_{2};
        dv=mydiff(v,wrt,indx_name);
        this=du+dv;
    case 'uplus'
        this=du;
    case 'minus'
        v=args_{2};
        dv=mydiff(v,wrt,indx_name);
        this=du-dv;
    case 'uminus'
        this=-du;
    case {'mtimes','times'}
        v=args_{2};
        dv=mydiff(v,wrt,indx_name);
        upv=du*v;
        vpu=dv*u;
        this=upv+vpu;
    case {'power','mpower'}
        v=args_{2};
        this=v*du*u^(v-1);
        if ~isnumeric(v.name)
            dv=mydiff(v,wrt,indx_name);
            this=this+dv*log(u)*u^v;
        end
    case {'rdivide','mrdivide'}
        v=args_{2};
        dv=mydiff(v,wrt,indx_name);
        this=(du*v-dv*u)/v^2; % only when v~=0
    case {'ldivide','mldivide'}
        v=args_{2};
        dv=mydiff(v,wrt,indx_name);
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
    otherwise
        error([myobj.name,' is unknown type of operator'])
end

end
