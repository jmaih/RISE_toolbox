classdef aplanar
    % aplanar automatic "planar" differentiation
    %
    % - [abs](aplanar/abs)
    % - [acos](aplanar/acos)
    % - [acosh](aplanar/acosh)
    % - [and](aplanar/and)
    % - [aplanar](aplanar/aplanar)
    % - [asin](aplanar/asin)
    % - [asinh](aplanar/asinh)
    % - [atan](aplanar/atan)
    % - [atanh](aplanar/atanh)
    % - [cos](aplanar/cos)
    % - [cosh](aplanar/cosh)
    % - [cot](aplanar/cot)
    % - [coth](aplanar/coth)
    % - [diff](aplanar/diff)
    % - [eq](aplanar/eq)
    % - [erf](aplanar/erf)
    % - [exp](aplanar/exp)
    % - [ge](aplanar/ge)
    % - [gt](aplanar/gt)
    % - [if_elseif](aplanar/if_elseif)
    % - [if_then_else](aplanar/if_then_else)
    % - [le](aplanar/le)
    % - [log](aplanar/log)
    % - [log10](aplanar/log10)
    % - [lt](aplanar/lt)
    % - [max](aplanar/max)
    % - [min](aplanar/min)
    % - [minus](aplanar/minus)
    % - [mpower](aplanar/mpower)
    % - [mrdivide](aplanar/mrdivide)
    % - [mtimes](aplanar/mtimes)
    % - [ne](aplanar/ne)
    % - [normcdf](aplanar/normcdf)
    % - [normpdf](aplanar/normpdf)
    % - [or](aplanar/or)
    % - [plus](aplanar/plus)
    % - [power](aplanar/power)
    % - [rdivide](aplanar/rdivide)
    % - [sign](aplanar/sign)
    % - [sin](aplanar/sin)
    % - [sinh](aplanar/sinh)
    % - [sqrt](aplanar/sqrt)
    % - [tan](aplanar/tan)
    % - [tanh](aplanar/tanh)
    % - [times](aplanar/times)
    % - [uminus](aplanar/uminus)
    % - [uplus](aplanar/uplus)
    properties
        x
        dx
        dxx
        dxxx
        dxxxx
        dxxxxx
        order
        nvars
    end
    methods
        % constructor
        %------------
        function obj=aplanar(v,nvars,vpos,order_)
            if nargin
                nv=numel(v);
                if nv~=numel(vpos)
                    error('number of positions should be the same as the number of locations')
                end
                if order_<=0||(order_~=floor(order_))
                    error('order of differentiation should be a strictly positive integer')
                end
                obj.order=order_;
                obj.nvars=nvars;
                further=repmat({nvars},1,order_);
                xxx=repmat('x',1,order_);
                for io=1:order_
                    dd=['d',xxx(1:io)];
                    ff=further(1:io);
                    obj.(dd)=zeros(1,ff{:});
                end
                if order_>1
                    obj.dxx=zeros(1,nvars,nvars);
                    if order_>2
                        obj.dxxx=zeros(1,nvars,nvars,nvars);
                    end
                end
                if ~isnan(vpos(1))
                    obj.dx(vpos(1))=1;
                end
                obj.x=v(1);
                if nv>1
                    tmp=obj;
                    for ii=2:nv
                        obj.x=v(ii);
                        obj.dx=zeros(1,nvars);
                        if ~isnan(vpos(ii))
                            obj.dx(vpos(ii))=1;
                        end
                        tmp(ii,1)=obj;
                    end
                    obj=tmp;
                end
            end
        end
        % unary
        %------
        function obj = abs(g)
            obj=apply_to_all(g,@times,sign(g.x));
        end
        function obj = acos(g)
            obj=0.5*pi-asin(g);
        end
        function obj = acosh(g)
            obj=2*log(sqrt(.5*(g+1))+sqrt(.5*(g-1)));
        end
        function obj = asin(g)
            obj=-1i*log(1i*g+sqrt(1-g^2));
        end
        function obj = asinh(g)
            obj=log(g+sqrt(1+g^2));
        end
        function obj = atan(g)
            obj=-1i*log((1+1i*g)*sqrt(1/(1+g^2)));
        end
        function obj = atanh(g)
            obj=.5*(log(1+g)-log(1-g));
        end
        function obj = cos(g)
            obj=g;
            sx=sin(g.x);
            cx=cos(g.x);
            obj.x=cx;
            obj.dx=-g.dx*sx;
            if obj.order>1
                gi_gj=aplanar.kron_gi_gj(g.dx,g.nvars);
                obj.dxx=-sx*g.dxx-cx*gi_gj;
                if obj.order>2
                    [sumperm_gijk,kron_gi_gj_gk]=aplanar.g_third(g.dx,g.dxx,g.nvars);
                    obj.dxxx=-(...
                        sx*(g.dxxx-kron_gi_gj_gk)+...
                        cx*sumperm_gijk...
                        );
                end
            end
        end
        function obj = cosh(g)
            obj=.5*(exp(g)+exp(-g));
        end
        function obj = cot(g)
            obj=1/tan(g);
        end
        function obj = coth(g)
            obj = cosh(g)/sinh(g);
        end
        function obj = erf(g)
            obj=g;
            obj.x=erf(g.x);
            obj.dx=2/sqrt(pi)*exp(-g.x^2)*g.dx;
            if obj.order>1
                gi_gj=aplanar.kron_gi_gj(g.dx,g.nvars);
                obj.dxx=2/sqrt(pi)*(g.dxx-2*gi_gj)*exp(-g.x^2);
                if obj.order>2
                    [sumperm_gijk,kron_gi_gj_gk]=aplanar.g_third(g.dx,g.dxx,g.nvars);
                    obj.dxxx=2/sqrt(pi)*(g.dxxx+4*kron_gi_gj_gk-2*sumperm_gijk)*exp(-g.x^2);
                end
            end
        end
        function obj = exp(g)
            obj=g;
            obj.x=exp(g.x);
            obj.dx=g.dx*obj.x;
            if obj.order>1
                gi_gj=aplanar.kron_gi_gj(g.dx,g.nvars);
                obj.dxx=(g.dxx+gi_gj)*obj.x;
                if obj.order>2
                    [sumperm_gijk,kron_gi_gj_gk]=aplanar.g_third(g.dx,g.dxx,g.nvars);
                    obj.dxxx=(g.dxxx+sumperm_gijk+kron_gi_gj_gk)*obj.x;
                end
            end
        end
        function obj = log(g)
            obj=g;
            obj.x=log(g.x);
            obj.dx=g.dx/g.x;
            if obj.order>1
                gi_gj=aplanar.kron_gi_gj(g.dx,g.nvars);
                obj.dxx=1/g.x*g.dxx-1/g.x^2*gi_gj;
                if obj.order>2
                    [sumperm_gijk,kron_gi_gj_gk]=aplanar.g_third(g.dx,g.dxx,g.nvars);
                    obj.dxxx=1/g.x*g.dxxx-1/g.x^2*sumperm_gijk+2/g.x^3*kron_gi_gj_gk;
                end
            end
        end
        function obj = log10(g)
            obj=log(g)/log(10);
        end
        function obj = sign(g)
            obj=zero_derivatives(@sign,g);
        end
        function obj = sin(g)
            obj=g;
            sx=sin(g.x);
            cx=cos(g.x);
            obj.x=sx;
            obj.dx=g.dx*cx;
            if obj.order>1
                gi_gj=aplanar.kron_gi_gj(g.dx,g.nvars);
                obj.dxx=cx*g.dxx-sx*gi_gj;
                if obj.order>2
                    [sumperm_gijk,kron_gi_gj_gk]=aplanar.g_third(g.dx,g.dxx,g.nvars);
                    obj.dxxx=cx*g.dxxx-sx*sumperm_gijk-cx*kron_gi_gj_gk;
                end
            end
        end
        function obj = sinh(g)
            obj=.5*(exp(g)-exp(-g));
        end
        function obj = sqrt(g)
            obj=g^0.5;
        end
        function obj = tan(g)
            obj = sin(g)/cos(g);
        end
        function obj = tanh(g)
            obj = sinh(g)/cosh(g);
        end
        function obj = uminus(g)
            obj=g;
            index=aplanar.indexes(obj.order);
            obj.x=-g.x;
            for ii=1:numel(index)
                ff=index{ii};
                obj.(ff)=-g.(ff);
            end
        end
        function obj = uplus(g)
            obj=g;
        end
        
        % binary
        %-------
        function obj = and(g,h)
            obj=zero_derivatives(@and,g,h);
        end
        function obj = eq(g,h)
            obj=zero_derivatives(@eq,g,h);
        end
        function obj = ge(g,h)
            obj=zero_derivatives(@ge,g,h);
        end
        function obj = gt(g,h)
            obj=zero_derivatives(@gt,g,h);
        end
        function obj = le(g,h)
            obj=zero_derivatives(@le,g,h);
        end
        function obj = lt(g,h)
            obj=zero_derivatives(@lt,g,h);
        end
        function obj = max(g,h)
            g_dbl=isa(g,'double');
            h_dbl=isa(h,'double');
            if g_dbl
                obj=h;
                if g>=obj.x
                    obj.x=g;
                    obj=zero_derivatives(@(x)x,obj);
                end
            elseif h_dbl
                obj=g;
                if h>=obj.x
                    obj.x=h;
                    obj=zero_derivatives(@(x)x,obj);
                end
            else
                if h.x>=g.x
                    obj=h;
                else
                    obj=g;
                end
            end
        end
        function obj = min(g,h)
            obj = -max(-g,-h);
        end
        function obj = minus(g,h)
            g_dbl=isa(g,'double');
            h_dbl=isa(h,'double');
            if g_dbl
                obj=-h;
                obj.x=obj.x+g;
            elseif h_dbl
                obj=g;
                obj.x=g.x-h;
            else
                obj=g;
                index=aplanar.indexes(obj.order);
                obj.x=g.x-h.x;
                for ii=1:numel(index)
                    ff=index{ii};
                    obj.(ff)=g.(ff)-h.(ff);
                end
            end
        end
        function obj = mpower(g,h)
            bad_flag=false;
            g_dbl=isa(g,'double');
            if ~g_dbl
                bad_flag=g.x<=0||~isreal(g.x);
                if g.x==0
                    g.x=eps*1i;
                end
            elseif g_dbl && g==0
                g=eps*1i;
                bad_flag=true;
            end
            obj=exp(h*log(g));
            if bad_flag
                myreal();
            end
            function myreal()
                if isa(obj,'double')
                    obj=real(obj);
                else
                    obj.x=real(obj.x);
                    xxx=repmat('x',1,obj.order);
                    for io=1:obj.order
                        dv=['d',xxx(1:io)];
                        obj.(dv)=real(obj.(dv));
                    end
                end
            end
        end
        function obj = mrdivide(g,h)
            obj=g*h^(-1);
        end
        function obj = mtimes(g,h)
            g_dbl=isa(g,'double');
            h_dbl=isa(h,'double');
            if g_dbl
                % if g==0,obj=0;
                obj=apply_to_all(h,@times,g);
            elseif h_dbl
                % if h==0,obj=0;
                obj=apply_to_all(g,@times,h);
            else
                obj=g;
                obj.x=g.x*h.x;
                obj.dx=g.x*h.dx+h.x*g.dx;
                if obj.order>1
                    obj.dxx=g.dxx*h.x+g.x*h.dxx;
                    for jj=1:obj.nvars
                        obj.dxx(1,:,jj)=obj.dxx(1,:,jj)+g.dx*h.dx(1,jj)+h.dx*g.dx(1,jj);
                    end
                    if obj.order>2
                        obj.dxxx=g.dxxx*h.x+g.x*h.dxxx;
                        for kk=1:obj.nvars
                            for jj=1:obj.nvars
                                obj.dxxx(1,:,jj,kk)=obj.dxxx(1,:,jj,kk)+...
                                    g.dxx(1,:,jj)*h.dx(1,kk)+...
                                    g.dxx(1,:,kk)*h.dx(1,jj)+...
                                    g.dx*h.dxx(1,jj,kk)+...
                                    h.dx*g.dxx(1,jj,kk)+...
                                    g.dx(1,jj)*h.dxx(1,:,kk)+...
                                    g.dx(1,kk)*h.dxx(1,:,jj);
                            end
                        end
                    end
                end
            end
        end
        function obj = ne(g,h)
            obj=zero_derivatives(@ne,g,h);
        end
        function obj = or(g,h)
            obj=zero_derivatives(@or,g,h);
        end
        function obj = plus(g,h)
            g_dbl=isa(g,'double');
            h_dbl=isa(h,'double');
            if g_dbl
                obj=h;
                obj.x=h.x+g;
            elseif h_dbl
                obj=g;
                obj.x=g.x+h;
            else
                obj=g;
                index=aplanar.indexes(obj.order);
                obj.x=g.x+h.x;
                for ii=1:numel(index)
                    ff=index{ii};
                    obj.(ff)=g.(ff)+h.(ff);
                end
            end
        end
        function obj = power(g,h)
            obj=mpower(g,h);
        end
        function obj = rdivide(g,h)
            obj=mrdivide(g,h);
        end
        function obj = times(g,h)
            obj=mtimes(g,h);
        end
        
        % n-ary
        %------
        function obj = if_elseif(varargin)
            % just pick the first element corresponding to the location
            % that evaluates to true
            nargs=length(varargin);
            if rem(nargs,2)
                error('the number of arguments must be even')
            end
            done=false;
            iarg=1;
            while ~done && iarg<nargs
                check=varargin{iarg};
                if ~isa(check,'double')
                    check=check.x;
                end
                if check
                    obj=varargin{iarg+1};
                    done=~done;
                end
                iarg=iarg+2;
            end
            if ~done
                error('could not find a valid statement in if_elseif')
            end
        end
        function obj = if_then_else(a,b,c)
            obj = if_elseif(a,b,1-a,c);
        end
        function obj = normcdf(x,mu,sig)
            if nargin<3
                sig=[];
                if nargin<2
                    mu=[];
                end
            end
            if isempty(sig)
                sig=1;
            end
            if isempty(mu)
                mu=0;
            end
            obj=0.5*(1+erf((x-mu)/(sig*sqrt(2))));
        end
        function obj = normpdf(x,mu,sig)
            if nargin<3
                sig=[];
                if nargin<2
                    mu=[];
                end
            end
            if isempty(sig)
                sig=1;
            end
            if isempty(mu)
                mu=0;
            end
            obj=exp(-.5*((x-mu)/sig)^2)/(sig*sqrt(2*pi));
        end
    end
    methods(Static,Access=private)
        function dxx=kron_gi_gj(gx,nvars)
            if nargin<2
                nvars=numel(gx);
            end
            %             tic
            dxx=zeros(1,nvars,nvars);
            for jj=1:nvars
                dxx(1,:,jj)=gx*gx(1,jj);
            end
            %             toc
            %             tic
            %             dxx2=bsxfun(@times,gx,permute(gx(1,:,ones(nvars,1)),[1,3,2]));
            %             toc
        end
        function [sumperm_gijk,kron_gi_gj_gk]=g_third(dx,dxx,nvars)
            if nargin<3
                nvars=numel(dx);
            end
            sumperm_gijk=zeros(1,nvars,nvars);
            kron_gi_gj_gk=sumperm_gijk;
            for kk=1:nvars
                for jj=1:nvars
                    sumperm_gijk(1,:,jj,kk)=dxx(1,:,kk)*dx(1,jj)+dx*dxx(1,jj,kk)+dxx(1,:,jj)*dx(1,kk);
                    kron_gi_gj_gk(1,:,jj,kk)=dx*dx(1,jj)*dx(1,kk);
                end
            end
        end
        function index=indexes(order)
            index=cell(1,order);
            index{1}='dx';
            for io=2:order
                index{io}=[index{io-1},'x'];
            end
        end
    end
    methods(Access=private)
        function obj=apply_to_all(obj,func,const)
            obj.x=func(obj.x,const);
            index=aplanar.indexes(obj.order);
            for ii=1:numel(index)
                ff=index{ii};
                obj.(ff)=func(obj.(ff),const);
            end
        end
        function obj=zero_derivatives(func,a,b)
            a_aplanar = isa(a,'aplanar');
            if nargin==3
                if a_aplanar
                    obj=a;
                    ax=a.x;
                    bx=b;
                else
                    obj=b;
                    bx=b.x;
                    ax=a;
                end
                args={ax,bx};
            elseif nargin==2
                obj=a;
                ax=a.x;
                args={ax};
            end
            obj.x=func(args{:});
            index=aplanar.indexes(obj.order);
            for ii=1:numel(index)
                obj.(index{ii})=0*obj.(index{ii});
            end
        end
    end
    methods(Static)
        function C=diff(func,active,inactive,order)
            n=numel(func);
            C=cell(1,order);
            nv=size(active,1);
            for ii=1:order
                C{ii}=zeros([n,nv*ones(1,ii)]);
            end
            silent=true;
            xxx=repmat('x',1,order);
            for ifunc=1:n
                % get info
                %---------
                objective=func{ifunc};
                arg_list=get_arg_list();
                % locate the variables for efficient differentiation
                %---------------------------------------------------
                pos_active=locate_variables(arg_list,active(:,1),silent);
                if isempty(inactive)
                    pos_inactive=nan*pos_active;
                else
                    pos_inactive=locate_variables(arg_list,inactive(:,1),silent);
                end
                nwrt=sum(~isnan(pos_active));
                nargs=numel(arg_list);
                vpos=0;
                x0=[];
                for iarg=1:nargs
                    if isnan(pos_active(iarg))
                        arg_list{iarg}=inactive{pos_inactive(iarg),2};
                    else
                        vpos=vpos+1;
                        if isempty(x0)
                            x0=aplanar(active{pos_active(iarg),2},nwrt,vpos,order);
                        else
                            x0.x=active{pos_active(iarg),2};
                            x0.dx(x0.dx~=0)=0;
                            x0.dx(vpos)=1;
                        end
                        arg_list{iarg}=x0;
                    end
                end
                % differentiate
                %--------------
                obj=objective(arg_list{:});
                % store the information
                %----------------------
                pos_active(isnan(pos_active))=[];
                if ~isempty(pos_active)
                    for io=1:order
                        further=repmat({pos_active},1,io);
                        C{io}(ifunc,further{:})=obj.(['d',xxx(1:io)]);
                    end
                end
            end
            function arg_list=get_arg_list()
                fstr=func2str(objective);
                right_par=find(fstr==')',1,'first');
                fstr=fstr(3:right_par-1);
                arg_list= regexp(fstr,',','split');
            end
        end
    end
end
