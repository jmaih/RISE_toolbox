function [xhandle,dxhandle,operation_number,operation_location,mycall] = differentiate(obj,wrt,wrt_index,fast)
% this function
if nargin<4
    fast=true;
end

switch class(wrt)
    case 'char'
        wrt=cellstr(wrt);
        wrt_order=1:numel(wrt);
    case 'cell'
        if isa(wrt{1},'sadiff')
            wrt=[wrt{:}];
            wrt_order=[wrt.order];
            wrt={wrt.func};
        else
             wrt_order=1:numel(wrt);
        end
    case 'sadiff'
        wrt_order=[wrt.order];
        wrt={wrt.func};
end
mycall=[];
F_dF_prefixes=[];
if fast
    mycall=sadiff.metastruct();
    F_dF_prefixes=mycall.prefix_list(2:3);
    indx_prefix=mycall.prefix_list{end};
end

nobj=numel(obj);
siz_obj=size(obj);
xhandle=cell(siz_obj);
dxhandle=cell(siz_obj);
operation_number=cell(siz_obj);
operation_location=cell(siz_obj);
for iobj=1:nobj
    if fast
        indx=mat2str(wrt_index{iobj});
        indx(isspace(indx))=',';
        indx_name=create_handle(iobj,indx_prefix);
        [~,~,mycall]=update_line(indx,indx_name,mycall);
    else
        indx_name=create_handle(iobj,'xx');
    end
    [xhandle{iobj},dxhandle{iobj},operation_number{iobj},operation_location{iobj},mycall] = mydiff(obj(iobj),fast,mycall,F_dF_prefixes);
end

mycall=sadiff.trim_metastruct(mycall,dxhandle);

    function [xhandle,dxhandle,operation_number,operation_location,mycall] = mydiff(obj_i,fast,mycall,F_dF_prefixes)
        
        write_definition=~isempty(mycall);
        if isa(obj_i.func,'double')
            % this is a constant
            val=sprintf('%0.10g',obj_i.func);
            der=@(z)'0';
            write_definition=false;
        elseif isempty(obj_i.args)
            % this is a variable name : use its order of
            % differentiation to define the derivative. if it has no
            % order then its derivative is 0
            val=obj_i.func;
            order_=wrt_order(strcmp(val,wrt));
            if isempty(order_)
                der=@(z)'0';
            else
                der=@(z)['bigi_(',sprintf('%0.10g',order_),',',indx_name,')'];
            end
            write_definition=false;
        else
            nargs=numel(obj_i.args);
            args_=struct('val',{},'der',{});
            for iarg=1:nargs
                % catch constants right here
                this_arg=obj_i.args{iarg};
                if isa(this_arg,'double')
                    this_arg=sadiff(this_arg);
                end
                [args_(iarg).val,args_(iarg).der,~,~,mycall]=mydiff(this_arg,fast,mycall,F_dF_prefixes);
            end
            % produce a final derivative out of those inputs
            [val,der]=get_derivative(obj_i.func,args_,nargs);
        end
        
        if write_definition
            mycall.operCount=mycall.operCount+1;
            [xhandle,dxhandle]=create_handle(mycall.operCount,F_dF_prefixes{:});
            [xhandle,xlineloc,mycall]=update_line(val,xhandle,mycall);
            [dxhandle,dxlineloc,mycall]=update_line(der(xhandle),dxhandle,mycall);
            operation_number=mycall.operCount;
            operation_location=[xlineloc,dxlineloc];
        else
            xhandle=val;
            dxhandle=der(val);
            operation_location=[];
            operation_number=[];
        end
        
        function [val,der]=get_derivative(func,args_,nargs)
            if nargs
                u=args_(1).val;
                du=args_(1).der;
                if nargs>1
                    v=args_(2).val;
                    dv=args_(2).der;
                    if nargs>2
                        w=args_(3).val;
                        %                     dw=args_(3).der;
                        if nargs>3
                            %                         x=args_(3).val;
                            %                         dx=args_(3).der;
                        end
                    end
                end
            end
            switch func
                % binary functions
                case {'^','power','mpower'}
                    val=sadiff.neat('^',u,v);
                    der21=sadiff.neat('*',v,du);
                    der22=sadiff.neat('-',v,'1');
                    der23=sadiff.neat('^',u,der22);
                    der2=sadiff.neat('*',der21,der23);
                    der11=sadiff.neat('*',dv,['log(',u,')']);
                    der1=@(z)sadiff.neat('*',der11,z);
                    der=@(z)sadiff.neat('+',der1(z),der2);
                case {'/','rdivide','mrdivide'}
                    val=sadiff.neat('/',u,v);
                    der1=sadiff.neat('*',du,v);
                    der2=sadiff.neat('*',dv,u);
                    der3=sadiff.neat('^',v,'2');
                    der=sadiff.neat('-',der1,der2);
                    der=@(z)sadiff.neat('/',der,der3);
                case {'\','ldivide','mldivide'}
                    [val,der]=get_derivative('/',args_(2:-1:1),nargs);
                case {'*','times','mtimes'}
                    val=sadiff.neat('*',u,v);
                    der1=sadiff.neat('*',du,v);
                    der2=sadiff.neat('*',dv,du);
                    der=@(z)sadiff.neat('+',der1,der2);
                case {'+','plus'}
                    val=sadiff.neat('+',u,v);
                    der=@(z)sadiff.neat('+',du,dv);
                case {'-','minus'}
                    val=sadiff.neat('-',u,v);
                    der=@(z)sadiff.neat('-',du,dv);
                    % unary functions
                case 'uminus'
                    val=sadiff.neat('-','0',u);
                    der=@(z)sadiff.neat('-','0',du);
                case 'uplus'
                    val=u;
                    der=@(z)du;
                case 'exp'
                    val=['exp(',u,')'];
                    der=@(z)sadiff.neat('*',du,z);
                case 'log'
                    val=['log(',u,')'];
                    der=@(z)sadiff.neat('/',du,u);
                case 'cos'
                    val=['cos(',u,')'];
                    der=sadiff.neat('*',du,['sin(',u,')']);
                    der=@(z)sadiff.neat('-','0',der);
                case 'acos'
                    val=['acos(',u,')'];
                    der0=sadiff.neat('^',u,'2');
                    der0=sadiff.neat('-','1',der0);
                    der0=sadiff.neat('/',du,['sqrt(',der0,')']);
                    der=@(z)sadiff.neat('-','0',der0);
                case 'cosh'
                    val=strcat('cosh(',u,')');
                    der=@(z)sadiff.neat('*',du,['sinh(',u,')']);
                case 'sin'
                    val=['sin(',u,')'];
                    der=@(z)sadiff.neat('*',du,['cos(',u,')']);
                case 'asin'
                    val=['asin(',u,')'];
                    der0=sadiff.neat('^',u,'2');
                    der0=sadiff.neat('-','1',der0);
                    der=@(z)sadiff.neat('/',du,['sqrt(',der0,')']);
                case 'sinh'
                    val=['sinh(',u,')'];
                    der=@(z)sadiff.neat('/',du,['cosh(',u,')']);
                case 'tan'
                    val=['tan(',u,')'];
                    der=sadiff.neat('^',['cos(',u,')'],'2');
                    der=@(z)sadiff.neat('/',du,der);
                case 'atan'
                    val=['atan(',u,')'];
                    der0=sadiff.neat('^',u,'2');
                    der0=sadiff.neat('+','1',der0);
                    der=@(z)sadiff.neat('/',du,['sqrt(',der0,')']);
                case 'tanh'
                    val=['tanh(',u,')'];
                    der0=sadiff.neat('^',['cosh(',u,')'],'2');
                    der=@(z)sadiff.neat('/',du,der0);
                case 'sqrt'
                    val=['sqrt(',u,')'];
                    der0=sadiff.neat('*','0.5',du);
                    der=@(z)sadiff.neat('/',der0,z);
                case 'abs'
                    val=['abs(',u,')'];
                    der=@(z)sadiff.neat('*',['sign(',u,')'],du);
                    % special functions
                case 'normalpdf'
                    if ~exist('v','var')
                        v='0';
                    end
                    if ~exist('w','var')
                        w='1';
                    end
                    sig=w;
                    mu=v;
                    val=['normpdf(',u,',',mu,',',sig,')'];
                    der0=sadiff.neat('-',mu,u);
                    der1=sadiff.neat('^',sig,'2');
                    der2=sadiff.neat('/',der0,der1);
                    der2=sadiff.neat('*',der2,du);
                    der = @(z)sadiff.neat('*',der2,z);
                case 'normalcdf'
                    if ~exist('v','var')
                        v='0';
                    end
                    if ~exist('w','var')
                        w='1';
                    end
                    sig=w;
                    mu=v;
                    val=['normcdf(',u,',',mu,',',sig,')'];
                    der0=sadiff.neat('*',du,['normpdf(',u,',',mu,',',sig,')']);
                    der = @(z)der0;
                otherwise
                    error(['method ',func,' is not implemented. please contact junior.maih@gmail.com'])
            end
            
        end
    end


end


