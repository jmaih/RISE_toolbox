classdef rise_nad
    % numerical automatic differentiation
    % see also sad, rise_sad
    % example objectives=@(x)[x(1)^2+x(2)*log(x(3));x(1)^2+x(2)*sin(x(3))];
    % vars=rise_nad(rand(3,1)); derivs=objectives(vars);
    % dx=vertcat(derivs.dx);
    
    % to do: unknown functions
    % jacobian
    % hessian
    properties
        x
        dx
    end
    methods
        function obj=rise_nad(y,dy)
            if nargin
                if isa(y,'rise_nad')
                    obj=y;
                    return
                end
                if nargin<2
                    dy=speye(length(y));
                end
                if numel(y)>1
                    obj=rise_nad.empty(0);
                    for iobj=1:numel(y)
                        obj(iobj,1).x=y(iobj);
                        obj(iobj,1).dx=dy(iobj,:);
                    end
                else
                    obj.x=y(1);
                    obj.dx=dy(1,:);
                end
            end
        end
        function obj=plus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            obj=rise_nad(u+v,du+dv);
        end
        function obj=uplus(u)
            [u,du]=get_props(u);
            obj=rise_nad(u,du);
        end
        function obj=minus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            obj=rise_nad(u-v,du-dv);
        end
        function obj=uminus(u)
            [u,du]=get_props(u);
            obj=rise_nad(-u,-du);
        end
        function obj=times(u,v)
            obj=mtimes(u,v);
        end
        function obj=mtimes(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            obj=rise_nad(u.*v,du*v+dv*u);
        end
        function obj=power(u,v)
            obj=mpower(u,v);
        end
        function obj=mpower(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=u.^v;
            if u==0
                n=max(length(du),length(dv));
                dxx=zeros(1,n);
            else
                dxx1=dv*log(u)*xx;
                dxx2=v.*du*u^(v-1);
                dxx=dxx1+dxx2;
            end
            obj=rise_nad(xx,dxx);
        end
        function obj=rdivide(u,v)
            obj=mrdivide(u,v);
        end
        function obj=mrdivide(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=u/v; 
            dxx=du.*v;
            dxx=dxx-dv*u;
            dxx=dxx/v^2;
            obj=rise_nad(xx,dxx);
        end
        function obj=ldivide(u,v)
            obj=mldivide(u,v);
        end
        function obj=mldivide(u,v)
            obj=mrdivide(v,u);
        end
        function obj=exp(u)
            [u,du]=get_props(u);
            xx=exp(u);
            dxx=du*xx;
            obj=rise_nad(xx,dxx);
        end
        function obj=log(u)
            [u,du]=get_props(u);
            xx=log(u);
            dxx=du./u;
            obj=rise_nad(xx,dxx);
        end
        function obj=log10(u)
            obj=log(u)/log(10);
        end
        function obj=cos(u)
            [u,du]=get_props(u);
            xx=cos(u);
            dxx=-du*sin(u);
            obj=rise_nad(xx,dxx);
        end
        function obj=acos(u)
            [u,du]=get_props(u);
            val=acos(u);
            der=sqrt(1-u.^2);
            der=-du/der;
            obj=rise_nad(val,der);
        end
        function obj=cosh(u)
            [u,du]=get_props(u);
            val=cosh(u);
            der=du*sinh(u);
            obj=rise_nad(val,der);
        end
        function obj=sin(u)
            [u,du]=get_props(u);
            xx=sin(u);
            dxx=du*cos(u);
            obj=rise_nad(xx,dxx);
        end
        function obj=asin(u)
            [u,du]=get_props(u);
            val=asin(u);
            der=sqrt(1-u^2);
            der=du./der;
            obj=rise_nad(val,der);
        end
        function obj=sinh(u)
            [u,du]=get_props(u);
            val=sinh(u);
            cu=cosh(u);
            der=du/cu;
            obj=rise_nad(val,der);
        end
        function obj=tan(u)
            [u,du]=get_props(u);
            val=tan(u);
            der=cos(u).^2;
            der=du/der;
            obj=rise_nad(val,der);
        end
        function obj=atan(u)
            [u,du]=get_props(u);
            val=atan(u);
            der=du/(1+u.^2);
            obj=rise_nad(val,der);
        end
        function obj=tanh(u)
            [u,du]=get_props(u);
            val=tanh(u);
            der=cosh(u)^2;
            der=du/der;
            obj=rise_nad(val,der);
        end
        function obj=min(u,v)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=u<v;
            val=min(u,v);
            der=uLv.*du+(1-uLv).*dv;
            obj = rise_nad(val,der);
        end
        function obj=max(u,v)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=u<v;
            val=max(u,v);
            der=uLv.*dv+(1-uLv).*du;
            obj = rise_nad(val,der);
        end
        function obj=sum(u,v)
            if nargin==1
                [u0,du0]=get_props(u(1));
                val=u0;
                der=du0;
                for ii=2:numel(u)
                    [u0,du0]=get_props(u(ii));
                    val=val+u0;
                    der=der+du0;
                end
                obj=rise_nad(val,der);
            else
                obj=plus(u,v);
            end
        end
        function obj=normpdf(u,mu,sig)
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            [u,du]=get_props(u);
            mu=get_props(mu);
            sig=get_props(sig);
            val=normpdf(u,mu,sig);
            der0=u-mu;
            der1=sig^2;
            der=-der0/der1;
            der=der.*du;
            der=der.*val;
            obj = rise_nad(val,der);
        end
        function obj=normcdf(u,mu,sig)
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            [u,du]=get_props(u);
            mu=get_props(mu);
            sig=get_props(sig);
            val=normcdf(u,mu,sig);
            der=du*normpdf(u,mu,sig);
            obj = rise_nad(val,der);
        end
        function obj=abs(u)
            [u,du]=get_props(u);
            val=abs(u);
            der=sign(u).*du;
            obj=rise_nad(val,der);
        end
        function obj=isreal(u)
            u=get_props(u);
            val=isreal(u);
            der=zeros(numel(u));
            obj=rise_nad(val,der);
        end
        function obj=sqrt(u)
            [u,du]=get_props(u);
            xx=sqrt(u);
            dxx=0.5.*du/xx;
            obj=rise_nad(xx,dxx);
        end
        function obj=norm(u)
            obj = sqrt(sum(u.^2));
        end
        function d=double(obj)
            d=vertcat(obj.dx);
        end
        %=====================
%         function a = subsasgn(a,s,b)
%             % SUBSASGN for rise_nad objects
%             if length(s)>1, error('Subscripted assignment too complicated'); end
%             if length(s.subs)==1, s.subs{end+1}=':'; end
%             
%             switch [class(a),class(b)]
%                 case 'doublerise_nad'
%                     a = subsasgn(a,s,double(b));
%                 case 'rise_naddouble'
%                     a.x  = subsasgn(a.x,s,b);
%                     a.dx = subsasgn(a.dx,s,0);
%                 case 'rise_nadrise_nad'
%                     a.x  = subsasgn(a.x,s,b.x);
%                     a.dx = subsasgn(full(a.dx),s,full(b.dx));
%                 otherwise
%                     error(['cannot assign ',class(b),' to ',class(a)]);
%             end
%         end
        %=====================
    end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
    end
end

function [u,du]=get_props(x)
du=0;
switch class(x)
    case 'rise_nad'
        u=x.x;
        du=x.dx;
    case 'double'
        u=x;
    otherwise
        error([mfilename,':: unsupported class ',class(x)])
end
end