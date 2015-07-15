classdef aplanar_
    % aplanar_ automatic "planar" differentiation
    %
    % - [abs](aplanar_/abs)
    % - [acos](aplanar_/acos)
    % - [acosh](aplanar_/acosh)
    % - [and](aplanar_/and)
    % - [aplanar_](aplanar_/aplanar_)
    % - [asin](aplanar_/asin)
    % - [asinh](aplanar_/asinh)
    % - [atan](aplanar_/atan)
    % - [atanh](aplanar_/atanh)
    % - [cos](aplanar_/cos)
    % - [cosh](aplanar_/cosh)
    % - [cot](aplanar_/cot)
    % - [coth](aplanar_/coth)
    % - [diff](aplanar_/diff)
    % - [eq](aplanar_/eq)
    % - [erf](aplanar_/erf)
    % - [exp](aplanar_/exp)
    % - [ge](aplanar_/ge)
    % - [gt](aplanar_/gt)
    % - [if_elseif](aplanar_/if_elseif)
    % - [if_then_else](aplanar_/if_then_else)
    % - [le](aplanar_/le)
    % - [log](aplanar_/log)
    % - [log10](aplanar_/log10)
    % - [lt](aplanar_/lt)
    % - [max](aplanar_/max)
    % - [min](aplanar_/min)
    % - [minus](aplanar_/minus)
    % - [mpower](aplanar_/mpower)
    % - [mrdivide](aplanar_/mrdivide)
    % - [mtimes](aplanar_/mtimes)
    % - [ne](aplanar_/ne)
    % - [normcdf](aplanar_/normcdf)
    % - [normpdf](aplanar_/normpdf)
    % - [or](aplanar_/or)
    % - [plus](aplanar_/plus)
    % - [power](aplanar_/power)
    % - [rdivide](aplanar_/rdivide)
    % - [sign](aplanar_/sign)
    % - [sin](aplanar_/sin)
    % - [sinh](aplanar_/sinh)
    % - [sqrt](aplanar_/sqrt)
    % - [tan](aplanar_/tan)
    % - [tanh](aplanar_/tanh)
    % - [times](aplanar_/times)
    % - [uminus](aplanar_/uminus)
    % - [uplus](aplanar_/uplus)
    properties
        x
        dx
        dxx
        dxxx
        dxxxx
        dxxxxx
        order
        nvars
        deriv_names
    end
    methods
        % constructor
        %------------
        function obj=aplanar_(v,nvars_,vpos,order_)
            % aplanar_ --  constructor for aplanar_ differentiation
            %
            % Syntax
            % -------
            % ::
            %
            %   obj=aplanar_(v,nvars,vpos,order_)
            %
            % Inputs
            % -------
            %
            % - **v** [vector]:
            %
            % - **nvars_** [integer]:
            %
            % - **vpos** [scalar|vector]:
            %
            % - **order_** [integer]:
            %
            % Outputs
            % --------
            %
            % - **obj** [aplanar_]: aplanar_ object with derivatives
            %
            % More About
            % ------------
            %
            % - it is desirable to further save some computations by
            % computing the shrinks only and then expanding entire maps as
            % in splanar...
            %
            % Examples
            % ---------
            %
            % See also:
            if nargin
                nv=numel(v);
                if nv~=numel(vpos)
                    error('number of positions should be the same as the number of locations')
                end
                if order_<=0||(order_~=floor(order_))
                    error('order of differentiation should be a strictly positive integer')
                end
                obj.order=order_;
                obj.nvars=nvars_;
                xxx=repmat('x',1,order_);
                dnames=cell(1,order_);
                for io=1:order_
                    dd=['d',xxx(1:io)];
                    obj.(dd)=zeros([1,nvars_*ones(1,io)]); % obj.(dd)=sparse(1,nchoosek(nvars_+io-1,io));
                    dnames{io}=dd;
                end
                obj.deriv_names=dnames;
                if ~isnan(vpos(1))
                    obj.dx(vpos(1))=1;
                end
                obj.x=v(1);
                if nv>1
                    tmp=obj;
                    for ii=2:nv
                        obj.x=v(ii);
                        obj.dx=sparse(1,nvars_);
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
            obj=sine_cosine(g,'cos');
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
            obj=reinitialize(g);
            f=erf(g.x);
            obj.x=f;
            oo=obj.order;
            for ii=1:g.nvars
                if ii==1
                    obj.dx=2/sqrt(pi)*exp(-g.x^2)*g.dx;
                    if oo==1
                        break
                    end
                end
                if oo>1 && obj.dx(ii)~=0
                    for jj=1:ii
                        do_second_order();
                        if oo>2 && obj.dxx(1,ii,jj)~=0
                            for kk=1:jj
                                do_third_order();
                                if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
                                    for ll=1:kk
                                        do_fourth_order();
                                        if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
                                            for mm=1:ll
                                                do_fifth_order();
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            function do_second_order()
                obj.dxx(1,ii,jj)=2*(g.dxx(1,ii,jj)*f-g.dx(jj)*obj.dx(ii));
            end
            function do_third_order()
                obj.dxxx(1,ii,jj,kk)=2*(g.dxxx(1,ii,jj,kk)*f...
                    +g.dxx(1,ii,jj)*obj.dx(kk)...
                    -g.dxx(1,jj,kk)*obj.dx(ii)...
                    -g.dx(jj)*g.dxx(1,ii,kk));
            end
            function do_fourth_order()
                obj.dxxxx(1,ii,jj,kk,ll)=2*(...
                    g.dxxx(1,ii,jj,kk,ll)*f...
                    +g.dxxx(1,ii,jj,kk)*obj.dx(ll)...
                    +g.dxxx(1,ii,jj,ll)*obj.dx(kk)...
                    +g.dxx(1,ii,jj)*obj.dxx(1,kk,ll)...
                    -g.dxxx(1,jj,kk,ll)*obj.dx(ii)...
                    -g.dxx(1,jj,kk)*obj.dxx(1,ii,ll)...
                    -g.dxx(1,jj,ll)*obj.dxx(1,ii,kk)...
                    -g.dx(jj)*obj.dxxx(1,ii,kk,ll)...
                    );
            end
            function do_fifth_order()
                obj.dxxxxx(1,ii,jj,kk,ll,mm)=2*(...
                    g.dxxxx(1,ii,jj,kk,ll,mm)*f...
                    +g.dxxxx(1,ii,jj,kk,ll)*obj.dx(mm)...
                    +g.dxxxx(1,ii,jj,kk,mm)*obj.dx(ll)...
                    +g.dxxx(1,ii,jj,kk)*obj.dxx(1,ll,mm)...
                    +g.dxxxx(1,ii,jj,ll,mm)*obj.dx(kk)...
                    +g.dxxx(1,ii,jj,ll)*obj.dxx(1,kk,mm)...
                    +g.dxxx(1,ii,jj,mm)*obj.dxx(1,kk,ll)...
                    +g.dxx(1,ii,jj)*obj.dxxx(1,kk,ll,mm)...
                    -g.dxxxx(1,jj,kk,ll,mm)*obj.dx(ii)...
                    -g.dxxx(1,jj,kk,ll)*obj.dxx(1,ii,mm)...
                    -g.dxxx(1,jj,kk,mm)*obj.dxx(1,ii,ll)...
                    -g.dxx(1,jj,kk)*obj.dxxx(1,ii,ll,mm)...
                    -g.dxxx(1,jj,ll,mm)*obj.dxx(1,ii,kk)...
                    -g.dxx(1,jj,ll)*obj.dxxx(1,ii,kk,mm)...
                    -g.dxx(1,jj,mm)*obj.dxxx(1,ii,kk,ll)...
                    -g.dx(jj)*obj.dxxxx(1,ii,kk,ll,mm)...
                    );
            end
        end
        function obj = exp(g)
            obj=reinitialize(g);
            f=exp(g.x);
            obj.x=f;
            oo=obj.order;
            for ii=1:g.nvars
                if ii==1
                    obj.dx=g.dx*f;
                    if oo==1
                        break
                    end
                end
                if oo>1 && obj.dx(ii)~=0
                    for jj=1:ii
                        do_second_order();
                        if oo>2 && obj.dxx(1,ii,jj)~=0
                            for kk=1:jj
                                do_third_order();
                                if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
                                    for ll=1:kk
                                        do_fourth_order();
                                        if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
                                            for mm=1:ll
                                                do_fifth_order();
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            function do_second_order()
                obj.dxx(1,ii,jj)=g.dxx(1,ii,jj)*f+g.dx(ii)*obj.dx(jj);
            end
            function do_third_order()
                obj.dxxx(1,ii,jj,kk)=g.dxxx(1,ii,jj,kk)*f+...
                    g.dxx(1,ii,jj)*obj.dx(kk)+...
                    g.dxx(1,ii,kk)*obj.dx(jj)+...
                    g.dx(ii)*obj.dxx(1,jj,kk);
            end
            function do_fourth_order()
                obj.dxxxx(1,ii,jj,kk,ll)=g.dxxxx(1,ii,jj,kk,ll)*f+...
                    g.dxxx(1,ii,jj,kk)*obj.dx(ll)+...
                    g.dxxx(1,ii,kk,ll)*obj.dx(jj)+...
                    g.dxx(1,ii,kk)*obj.dxx(1,jj,ll)+...
                    g.dxxx(1,jj,kk,ll)*obj.dx(ii)+...
                    g.dxx(1,jj,kk)*obj.dxx(1,ii,ll)+...
                    g.dxx(1,kk,ll)*obj.dxx(1,ii,jj)+...
                    g.dx(kk)*obj.dxxx(1,ii,jj,ll);
            end
            function do_fifth_order()
                obj.dxxxxx(1,ii,jj,kk,ll,mm)=g.dxxxxx(1,ii,jj,kk,ll,mm)*f+...
                    g.dxxxx(1,ii,jj,kk,ll)*obj.dx(mm)+...
                    g.dxxxx(1,ii,jj,kk,mm)*obj.dx(ll)+...
                    g.dxxxx(1,ii,kk,ll,mm)*obj.dx(jj)+...
                    g.dxxxx(1,jj,kk,ll,mm)*obj.dx(ii)+...
                    g.dxxx(1,ii,jj,kk)*obj.dxx(1,ll,mm)+...
                    g.dxxx(1,ii,kk,ll)*obj.dxx(1,jj,mm)+...
                    g.dxxx(1,ii,kk,mm)*obj.dxx(1,jj,ll)+...
                    g.dxxx(1,jj,kk,ll)*obj.dxx(1,ii,mm)+...
                    g.dxxx(1,jj,kk,mm)*obj.dxx(1,ii,ll)+...
                    g.dxxx(1,kk,ll,mm)*obj.dxx(1,ii,jj)+...
                    g.dxx(1,ii,kk)*obj.dxxx(1,jj,ll,mm)+...
                    g.dxx(1,jj,kk)*obj.dxxx(1,ii,ll,mm)+...
                    g.dxx(1,kk,ll)*obj.dxxx(1,ii,jj,mm)+...
                    g.dxx(1,kk,mm)*obj.dxxx(1,ii,jj,ll)+...
                    g.dx(kk)*obj.dxxxx(1,ii,jj,ll,mm);
            end
        end
        function obj = log(g)
            obj=reinitialize(g);
            obj.x=log(g.x);
            oo=obj.order;
            igx=1/g.x;
            for ii=1:g.nvars
                if ii==1
                    obj.dx=g.dx*igx;
                    if oo==1
                        break
                    end
                end
                if oo>1 && obj.dx(ii)~=0
                    for jj=1:ii
                        do_second_order();
                        if oo>2 && obj.dxx(1,ii,jj)~=0
                            for kk=1:jj
                                do_third_order();
                                if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
                                    for ll=1:kk
                                        do_fourth_order();
                                        if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
                                            for mm=1:ll
                                                do_fifth_order();
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            function do_second_order()
                obj.dxx(1,ii,jj)=igx*g.dxx(1,ii,jj)-obj.dx(ii)*obj.dx(jj);
            end
            function do_third_order()
                obj.dxxx(1,ii,jj,kk)=igx*(g.dxxx(1,ii,jj,kk)-g.dxx(1,ii,jj)*obj.dx(kk))...
                    -(...
                    obj.dxx(1,jj,kk)*obj.dx(ii)...
                    +obj.dxx(1,ii,kk)*obj.dx(jj)...
                    );
            end
            function do_fourth_order()
                obj.dxxxx(1,ii,jj,kk,ll)=igx*(...
                    g.dxxxx(1,ii,jj,kk,ll)...
                    -g.dxxx(1,ii,jj,kk)*obj.dx(ll)...
                    -(g.dxxx(1,ii,jj,ll)-g.dxx(1,ii,jj)*obj.dx(ll))*obj.dx(kk)...
                    -g.dxx(1,ii,jj)*obj.dxx(1,kk,ll)...
                    )-(...
                    obj.dxxx(1,jj,kk,ll)*obj.dx(ii)...
                    +obj.dxxx(1,ii,kk,ll)*obj.dx(jj)...
                    +obj.dxx(1,jj,kk)*obj.dxx(1,ii,ll)...
                    +obj.dxx(1,jj,ll)*obj.dxx(1,ii,kk));
            end
            function do_fifth_order()
                a1=g.dxxxxx(1,ii,jj,kk,ll,mm)-g.dxxxx(1,ii,jj,kk,ll)*obj.dx(mm);
                
                a2=(g.dxxxx(1,ii,jj,kk,mm)-g.dxxx(1,ii,jj,kk)*obj.dx(mm))*obj.dx(ll)+...
                    g.dxxx(1,ii,jj,kk)*obj.dxx(1,ll,mm);
                
                a3=(g.dxxxx(1,ii,jj,ll,mm)-g.dxxx(1,ii,jj,ll)*obj.dx(mm))*obj.dx(kk)+...
                    g.dxxx(1,ii,jj,ll)*obj.dxx(1,kk,mm);
                
                a4=(g.dxxx(1,ii,jj,mm)-g.dxx(1,ii,jj)*obj.dx(mm))*obj.dx(ll)*obj.dx(kk)+...
                    g.dxx(1,ii,jj)*obj.dxx(1,ll,mm)*obj.dx(kk)+...
                    g.dxx(1,ii,jj)*obj.dx(ll)*obj.dxx(1,kk,mm);
                
                a5=(g.dxxx(1,ii,jj,mm)-g.dxx(1,ii,jj)*obj.dx(mm))*obj.dxx(1,kk,ll)+...
                    g.dxx(1,ii,jj)*obj.dxxx(1,kk,ll,mm);
                
                a6=obj.dxxxx(1,ii,kk,ll,mm)*obj.dx(jj)+...
                    obj.dxxx(1,ii,kk,ll)*obj.dxx(1,jj,mm)+...
                    obj.dxxx(1,ii,kk,mm)*obj.dxx(1,jj,ll)+...
                    obj.dxx(1,ii,kk)*obj.dxxx(1,jj,ll,mm)+...
                    obj.dxxx(1,ii,ll,mm)*obj.dxx(1,jj,kk)+...
                    obj.dxx(1,ii,ll)*obj.dxxx(1,jj,kk,mm)+...
                    obj.dxx(1,ii,mm)*obj.dxxx(1,jj,kk,ll)+...
                    obj.dx(ii)*obj.dxxxx(1,jj,kk,ll,mm);
                
                obj.dxxxxx(1,ii,jj,kk,ll,mm)=igx*(a1-a2-a3+a4-a5)-a6;
            end
        end
        function obj = log10(g)
            obj=log(g)/log(10);
        end
        function obj = sign(g)
            obj=zero_derivatives(@sign,g);
        end
        function obj = sin(g)
            obj=sine_cosine(g,'sin');
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
            index=obj.deriv_names;
            obj.x=-g.x;
            for ii=1:obj.order
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
                index=obj.deriv_names;
                obj.x=g.x-h.x;
                for ii=1:obj.order
                    ff=index{ii};
                    obj.(ff)=g.(ff)-h.(ff);
                end
            end
        end
        function obj = mpower(g,h)
            g_dbl=isa(g,'double');
            if ~g_dbl
                bad_flag=g.x<=0||~isreal(g.x);
                item=g;
            elseif g_dbl && g==0
                bad_flag=true;
                item=h;
            end
            if bad_flag
                obj=reinitialize(item);
            else
                obj=exp(h*log(g));
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
                obj=reinitialize(g);
                obj.x=g.x*h.x;
                oo=obj.order;
                for ii=1:g.nvars
                    if ii==1
                        obj.dx=g.x*h.dx+h.x*g.dx;
                        if oo==1
                            break
                        end
                    end
                    if oo>1 && obj.dx(ii)~=0
                        for jj=1:ii
                            do_second_order();
                            if oo>2 && obj.dxx(1,ii,jj)~=0
                                for kk=1:jj
                                    do_third_order();
                                    if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
                                        for ll=1:kk
                                            do_fourth_order();
                                            if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
                                                for mm=1:ll
                                                    do_fifth_order();
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            function do_second_order()
                obj.dxx(1,ii,jj)=g.dxx(1,ii,jj)*h.x...
                    +g.dx(ii)*h.dx(jj)...
                    +g.dx(jj)*h.dx(ii)...
                    +g.x*h.dxx(1,ii,jj);
            end
            function do_third_order()
                obj.dxxx(1,ii,jj,kk)=g.dxxx(1,ii,jj,kk)*h.x...
                    +g.dxx(1,ii,jj)*h.dx(kk)...
                    +g.dxx(1,ii,kk)*h.dx(jj)...
                    +g.dxx(1,jj,kk)*h.dx(ii)...
                    +g.dx(ii)*h.dxx(1,jj,kk)...
                    +g.dx(jj)*h.dxx(1,ii,kk)...
                    +g.dx(kk)*h.dxx(1,ii,jj)...
                    +g.x*h.dxxx(1,ii,jj,kk);
            end
            function do_fourth_order()
                obj.dxxxx(1,ii,jj,kk,ll)=g.dxxxx(1,ii,jj,kk,ll)*h.x...
                    +g.dxxx(1,ii,jj,kk)*h.dx(ll)...
                    +g.dxxx(1,ii,jj,ll)*h.dx(kk)...
                    +g.dxxx(1,ii,kk,ll)*h.dx(jj)...
                    +g.dxxx(1,jj,kk,ll)*h.dx(ii)...
                    +g.dxx(1,ii,jj)*h.dxx(1,kk,ll)...
                    +g.dxx(1,ii,kk)*h.dxx(1,jj,ll)...
                    +g.dxx(1,ii,ll)*h.dxx(1,jj,kk)...
                    +g.dxx(1,jj,kk)*h.dxx(1,ii,ll)...
                    +g.dxx(1,jj,ll)*h.dxx(1,ii,kk)...
                    +g.dxx(1,kk,ll)*h.dxx(1,ii,jj)...
                    +g.dx(jj)*h.dxxx(1,ii,kk,ll)...
                    +g.dx(kk)*h.dxxx(1,ii,jj,ll)...
                    +g.dx(ll)*h.dxxx(1,ii,jj,kk)...
                    +g.dx(ii)*h.dxxx(1,jj,kk,ll)...
                    +g.x*h.dxxxx(1,ii,jj,kk,ll);
            end
            function do_fifth_order()
                obj.dxxxxx(1,ii,jj,kk,ll,mm)=g.dxxxxx(1,ii,jj,kk,ll,mm)*h.x...
                    +g.dxxxx(1,ii,jj,kk,ll)*h.dx(mm)...
                    +g.dxxxx(1,ii,jj,kk,mm)*h.dx(ll)...
                    +g.dxxxx(1,ii,kk,ll,mm)*h.dx(jj)...
                    +g.dxxxx(1,jj,kk,ll,mm)*h.dx(ii)...
                    +g.dxxxx(1,ii,jj,ll,mm)*h.dx(kk)...
                    +g.dxxx(1,ii,jj,kk)*h.dxx(1,ll,mm)...
                    +g.dxxx(1,ii,jj,mm)*h.dxx(1,kk,ll)...
                    +g.dxxx(1,ii,kk,ll)*h.dxx(1,jj,mm)...
                    +g.dxxx(1,ii,kk,mm)*h.dxx(1,jj,ll)...
                    +g.dxxx(1,ii,jj,ll)*h.dxx(1,kk,mm)...
                    +g.dxxx(1,ii,ll,mm)*h.dxx(1,jj,kk)...
                    +g.dxxx(1,jj,kk,mm)*h.dxx(1,ii,ll)...
                    +g.dxxx(1,jj,ll,mm)*h.dxx(1,ii,kk)...
                    +g.dxxx(1,kk,ll,mm)*h.dxx(1,ii,jj)...
                    +g.dxxx(1,jj,kk,ll)*h.dxx(1,ii,mm)...
                    +g.dxx(1,ii,jj)*h.dxxx(1,kk,ll,mm)...
                    +g.dxx(1,ii,kk)*h.dxxx(1,jj,ll,mm)...
                    +g.dxx(1,ii,ll)*h.dxxx(1,jj,kk,mm)...
                    +g.dxx(1,ii,mm)*h.dxxx(1,jj,kk,ll)...
                    +g.dxx(1,jj,kk)*h.dxxx(1,ii,ll,mm)...
                    +g.dxx(1,jj,ll)*h.dxxx(1,ii,kk,mm)...
                    +g.dxx(1,jj,mm)*h.dxxx(1,ii,kk,ll)...
                    +g.dxx(1,kk,ll)*h.dxxx(1,ii,jj,mm)...
                    +g.dxx(1,kk,mm)*h.dxxx(1,ii,jj,ll)...
                    +g.dxx(1,ll,mm)*h.dxxx(1,ii,jj,kk)...
                    +g.dx(ii)*h.dxxxx(1,jj,kk,ll,mm)...
                    +g.dx(jj)*h.dxxxx(1,ii,kk,ll,mm)...
                    +g.dx(kk)*h.dxxxx(1,ii,jj,ll,mm)...
                    +g.dx(ll)*h.dxxxx(1,ii,jj,kk,mm)...
                    +g.dx(mm)*h.dxxxx(1,ii,jj,kk,ll)...
                    +g.x*h.dxxxxx(1,ii,jj,kk,ll,mm);
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
                index=obj.deriv_names;
                obj.x=g.x+h.x;
                for ii=1:obj.order
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
    
    methods(Access=private)
        function obj=apply_to_all(obj,func,const)
            obj.x=func(obj.x,const);
            index=obj.deriv_names;
            for ii=1:obj.order
                ff=index{ii};
                obj.(ff)=func(obj.(ff),const);
            end
        end
        function obj=zero_derivatives(func,a,b)
            a_aplanar = isa(a,'aplanar_');
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
            index=obj.deriv_names;
            for ii=1:obj.order
                obj.(index{ii})(:)=0;
            end
        end
    end
    methods(Static)
        function C=diff(func,active,inactive,order)
            n=numel(func);
            nv=size(active,1);
            recorder=struct('iter',{},'ncells',{},'info',{});
            ncells=1000;
            for io=1:order
                recorder(io).iter=0;
                recorder(io).ncells=ncells;
                recorder(io).info=cell(ncells,1);
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
                            x0=aplanar_(active{pos_active(iarg),2},nwrt,vpos,order);
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
                    % important that the following be horizontal
                    %--------------------------------------------
                    pos_active=pos_active(:).';
                    footprint=(1:nwrt).';
                    for io=1:order
                        derivs=obj.(['d',xxx(1:io)]);
                        if io==1
                            good=find(derivs);
                            posisjon=pos_active(good);
                            record_derivatives([],derivs(good),posisjon);
                        else
                            col_pages_books_etc=fliplr(2:io+1);
                            derivs=permute(derivs,[1,col_pages_books_etc]);
                            derivs=reshape(derivs,1,[]);
                            
                            prototype=cell2mat(utils.gridfuncs.mypermutation(1:io));
                            contract=utils.kronecker.shrink_expand(nwrt,io);
                            
%                             raw=utils.gridfuncs.mygrid(nwrt*ones(1,io));
%                             disp('max(max(raw(contract,:)-footprint))')
%                             max(max(raw(contract,:)-footprint))
%                             old_derivs=derivs;
                            
                            derivs=derivs(contract);
                            good=find(derivs);
                            derivs=derivs(good);
                            myprints=pos_active(footprint(good,:));
                            for id=1:numel(derivs)
                                record_derivatives(myprints(id,:),derivs(id));
                            end
                        end
                        if io<order
                            footprint=utils.gridfuncs.derivatives_grid(footprint);
                        end
                    end
                end
            end
            % now create output matrices
            %----------------------------
            C=cell(1,order);
            for io=1:order
                info=cell2mat(recorder(io).info(1:recorder(io).iter));
                C{io}=sparse(info(:,1),info(:,2),info(:,3),n,nv^io);
            end
            
            function pos=record_derivatives(array,d,pos)
                if nargin<3
                    wrt_wise=false;
                    pos0=array(prototype);
                    pos=utils.gridfuncs.locate_permutation(pos0,nv,wrt_wise);
                    pos=unique(pos);
                end
                iter=recorder(io).iter+1;
                recorder(io).iter=recorder(io).iter+1;
                if iter>=recorder(io).ncells
                    recorder(io).ncells=recorder(io).ncells+ncells;
                    recorder(io).info{end+ncells}={};
                end
                npos=numel(pos);
                if isscalar(d)
                    d=d(ones(npos,1));
                end
                recorder(io).info{iter}=[ifunc(ones(npos,1)),pos(:),d(:)];
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

function obj=sine_cosine(g,type)
switch type
    case 'cos'
        f=cos(g.x);
        z=-sin(g.x);
    case 'sin'
        f=sin(g.x);
        z=cos(g.x);
    otherwise
        error(['unknown type ',type])
end
obj=reinitialize(g);

obj.x=f;
oo=obj.order;
for ii=1:obj.nvars
    if ii==1
        obj.dx=g.dx*z;
        if oo==1
            break
        end
    end
    if oo>1 && obj.dx(ii)~=0
        for jj=1:ii
            do_second_order();
            if oo>2 && obj.dxx(1,ii,jj)~=0
                for kk=1:jj
                    do_third_order();
                    if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
                        for ll=1:kk
                            do_fourth_order();
                            if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
                                for mm=1:ll
                                    do_fifth_order();
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

    function do_second_order()
        obj.dxx(1,ii,jj)=z*g.dxx(1,ii,jj)-g.dx(ii)*g.dx(jj)*f;
    end
    function do_third_order()
        obj.dxxx(1,ii,jj,kk)=-(g.dx(kk)*g.dxx(1,ii,jj)+...
            g.dxx(1,jj,kk)*g.dx(ii)+...
            g.dx(jj)*g.dxx(1,ii,kk))*f...
            +z*g.dxxx(1,ii,jj,kk)...
            -g.dx(ii)*g.dx(jj)*obj.dx(kk);
    end
    function do_fourth_order()
        obj.dxxxx(1,ii,jj,kk,ll)=-(...
            g.dxx(1,kk,ll)*g.dxx(1,ii,jj)+...
            g.dx(kk)*g.dxxx(1,ii,jj,ll)+...
            g.dxxx(1,jj,kk,ll)*g.dx(ii)+...
            g.dxx(1,jj,kk)*g.dxx(1,ii,ll)+...
            g.dxx(1,jj,ll)*g.dxx(1,ii,kk)+...
            g.dx(jj)*g.dxxx(1,ii,kk,ll)+...
            g.dx(ll)*g.dxxx(1,ii,jj,kk)...
            )*f...
            -(...
            g.dx(kk)*g.dxx(1,ii,jj)+...
            g.dxx(1,jj,kk)*g.dx(ii)+...
            g.dx(jj)*g.dxx(1,ii,kk)...
            )*obj.dx(ll)...
            -(...
            g.dxx(1,jj,ll)*g.dx(ii)+...
            g.dx(jj)*g.dxx(1,ii,ll)...
            )*obj.dx(kk)...
            -g.dx(jj)*g.dx(ii)*obj.dxx(1,kk,ll)...
            +z*g.dxxxx(1,ii,jj,kk,ll);
    end
    function do_fifth_order()
        obj.dxxxxx(1,ii,jj,kk,ll,mm)=...
            z*g.dxxxxx(1,ii,jj,kk,ll,mm)-(...
            (...
            g.dx(mm)*g.dxxxx(1,ii,jj,kk,ll)+...
            g.dxxx(1,kk,ll,mm)*g.dxx(1,ii,jj)+g.dxx(1,kk,ll)*g.dxxx(1,ii,jj,mm)+...
            g.dxx(1,kk,mm)*g.dxxx(1,ii,jj,ll)+g.dx(kk)*g.dxxxx(1,ii,jj,ll,mm)+...
            g.dxxxx(1,jj,kk,ll,mm)*g.dx(ii)+g.dxxx(1,jj,kk,ll)*g.dxx(1,ii,mm)+...
            g.dxxx(1,jj,kk,mm)*g.dxx(1,ii,ll)+g.dxx(1,jj,kk)*g.dxxx(1,ii,ll,mm)+...
            g.dxxx(1,jj,ll,mm)*g.dxx(1,ii,kk)+g.dxx(1,jj,ll)*g.dxxx(1,ii,kk,mm)+...
            g.dxx(1,jj,mm)*g.dxxx(1,ii,kk,ll)+g.dx(jj)*g.dxxxx(1,ii,kk,ll,mm)+...
            g.dxx(1,ll,mm)*g.dxxx(1,ii,jj,kk)+g.dx(ll)*g.dxxxx(1,ii,jj,kk,mm)...
            )*f...
            +(g.dxx(1,kk,ll)*g.dxx(1,ii,jj)+...
            g.dx(kk)*g.dxxx(1,ii,jj,ll)+...
            g.dxxx(1,jj,kk,ll)*g.dx(ii)+...
            g.dxx(1,jj,kk)*g.dxx(1,ii,ll)+...
            g.dxx(1,jj,ll)*g.dxx(1,ii,kk)+...
            g.dx(jj)*g.dxxx(1,ii,kk,ll)+...
            g.dx(ll)*g.dxxx(1,ii,jj,kk))*obj.dx(mm)...
            +(...
            g.dxx(1,kk,mm)*g.dxx(1,ii,jj)+g.dx(kk)*g.dxxx(1,ii,jj,mm)+...
            g.dxxx(1,jj,kk,mm)*g.dx(ii)+g.dxx(1,jj,kk)*g.dxx(1,ii,mm)+...
            g.dxx(1,jj,mm)*g.dxx(1,ii,kk)+g.dx(jj)*g.dxxx(1,ii,kk,mm)...
            )*obj.dx(ll)...
            +(...
            g.dx(kk)*g.dxx(1,ii,jj)+...
            g.dxx(1,jj,kk)*g.dx(ii)+...
            g.dx(jj)*g.dxx(1,ii,kk)...
            )*obj.dxx(1,ll,mm)...
            +(...
            g.dxxx(1,jj,ll,mm)*g.dx(ii)+g.dxx(1,jj,ll)*g.dxx(1,ii,mm)+...
            g.dxx(1,jj,mm)*g.dxx(1,ii,ll)+g.dx(jj)*g.dxxx(1,ii,ll,mm)...
            )*obj.dx(kk)...
            +(...
            g.dxx(1,jj,ll)*g.dx(ii)+...
            g.dx(jj)*g.dxx(1,ii,ll)...
            )*obj.dxx(1,kk,mm)...
            +g.dxx(1,jj,mm)*g.dx(ii)*obj.dxx(1,kk,ll)...
            +g.dx(jj)*g.dxx(1,ii,mm)*obj.dxx(1,kk,ll)...
            +g.dx(jj)*g.dx(ii)*obj.dxxx(1,kk,ll,mm)...
            );
    end
end

function obj=reinitialize(g)
obj=g;
obj.x=0;
oo_=g.order;
obj.dx(:)=0;
if oo_>1
    obj.dxx(:)=0;
    if oo_>2
        obj.dxxx(:)=0;
        if oo_>3
            obj.dxxxx(:)=0;
            if oo_>4
                obj.dxxxxx(:)=0;
            end
        end
    end
end
end

%{
clc
fid = fopen('aplanar.m');A=char(fread(fid)');fclose(fid);

A=regexprep(A,'\(1,(\w{2},)+(\w{2})\)','(dl([$1$2]))');
A=regexprep(A,'\(1,(\w{2})\)','($1)');
A
%%
test='g.dxxx(1,kk,ll,mm)*g.dxx(1,ii,jj)+g.dxx(1,kk,ll)*g.dxxx(1,ii,jj,mm)+...';
regexprep(test,'\(1,(\w{2},)+(\w{2})\)','(dl([$1$2]))')

%-----------------

clear all
loc='/Users/juniormaih/Documents/RISE_toolbox/m/differentiation/@aplanar_';
fid=fopen([loc,filesep,'aplanar_.m'],'r');
c=char(fread(fid).');
fclose(fid);
%%
clc
test=c;%'g.dxx(dl([kk,mm]))*g.dxx(dl([ii,jj]))+g.dx(kk)*g.dxxx(dl([ii,jj,mm]))';
A=regexprep(test,'dl\(\[((?>\w{2},)+)(\w{2})\]\)','1,$1$2')

%}