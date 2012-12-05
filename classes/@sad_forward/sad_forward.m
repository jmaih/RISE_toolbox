classdef sad_forward %< handle
    % to do: unknown functions
    properties
        x
        dx
    end
    methods
        function obj=sad_forward(y,dy)
            if nargin
                if isa(y,'sad_forward')
                    obj=y;
                    return
                end
                if nargin<2
                    dy='0';
                end
                if ~ischar(y)||~ischar(dy)
                    error([mfilename,':: inputs must be char'])
                end
                obj.x=y;
                obj.dx=dy;
            end
        end
        function obj=plus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=rise_algebra_cat('+',u,v);   
            dxx=rise_algebra_cat('+',du,dv);
            obj=sad_forward(xx,dxx);
        end
        function obj=uplus(u)
            [u,du]=get_props(u);
            xx=u;
            dxx=du;
            obj=sad_forward(xx,dxx);
        end
        function obj=minus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=rise_algebra_cat('-',u,v);    
            dxx=rise_algebra_cat('-',du,dv); 
            obj=sad_forward(xx,dxx);
        end
        function obj=uminus(u)
            [u,du]=get_props(u);
            xx=rise_algebra_cat('-',u);  
            dxx=rise_algebra_cat('-',du);
            obj=sad_forward(xx,dxx);
        end
        function obj=times(u,v)
            obj=mtimes(u,v);
        end
        function obj=mtimes(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=rise_algebra_cat('*',u,v); 
            dxx1=rise_algebra_cat('*',du,v);
            dxx2=rise_algebra_cat('*',dv,u);
            dxx=rise_algebra_cat('+',dxx1,dxx2);
            obj=sad_forward(xx,dxx);
        end
        function obj=power(u,v)
            obj=mpower(u,v);
        end
        function obj=mpower(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=rise_algebra_cat('^',u,v);
            dxx1=rise_algebra_cat('*',rise_algebra_cat('*',dv,['log(',u,')']),xx);
            dxx2=rise_algebra_cat('*',rise_algebra_cat('*',v,du),...
                rise_algebra_cat('^',u,rise_algebra_cat('-',v,'1')));
            dxx=rise_algebra_cat(dxx1,dxx2,'+');
            obj=sad_forward(xx,dxx);
        end
        function obj=rdivide(u,v)
            obj=mrdivide(u,v);
        end
        function obj=mrdivide(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=rise_algebra_cat('/',u,v);
            dxx=rise_algebra_cat('*',du,v);
            dxx=rise_algebra_cat('-',dxx,rise_algebra_cat('*',dv,u));
            dxx=rise_algebra_cat('/',dxx,rise_algebra_cat('^',v,'2'));
            obj=sad_forward(xx,dxx);
        end
        function obj=ldivide(u,v)
            obj=mldivide(u,v);
        end
        function obj=mldivide(u,v)
            obj=mrdivide(v,u);
        end
        function obj=exp(u)
            [u,du]=get_props(u);
            xx=['exp(',u,')'];
            dxx=rise_algebra_cat('*',du,xx);
            obj=sad_forward(xx,dxx);
        end
        function obj=log(u)
            [u,du]=get_props(u);
            xx=['log(',u,')'];
            dxx=rise_algebra_cat('/',du,u);
            obj=sad_forward(xx,dxx);
        end
        function obj=log10(u)
            obj=log(u)/log(10);
        end
        function obj=cos(u)
            [u,du]=get_props(u);
            xx=['cos(',u,')'];
            dxx=rise_algebra_cat('*',...
                rise_algebra_cat('-','0',du),['sin(',u,')']);
            obj=sad_forward(xx,dxx);
        end
        function obj=acos(u)
            [u,du]=get_props(u);
            val=['acos(',u,')'];
            der=rise_algebra_cat('^',u,'2');
            der=rise_algebra_cat('-','1',der);
            der=rise_algebra_cat('/',rise_algebra_cat('*','-1',du),...
                ['sqrt(',der,')']);
            obj=sad_forward(val,der);
        end
        function obj=cosh(u)
            [u,du]=get_props(u);
            val=['cosh(',u,')'];
            der=rise_algebra_cat('*',du,['sinh(',u,')']);
            obj=sad_forward(val,der);
        end
        function obj=sin(u)
            [u,du]=get_props(u);
            xx=['sin(',u,')'];
            dxx=rise_algebra_cat('*',du,['cos(',u,')']);
            obj=sad_forward(xx,dxx);
        end
        function obj=asin(u)
            [u,du]=get_props(u);
            val=['asin(',u,')'];
            der=rise_algebra_cat('^',u,'2');
            der=rise_algebra_cat('*','-1',der);
            der=rise_algebra_cat('/',du,['sqrt(',der,')']);
            obj=sad_forward(val,der);
        end
        function obj=sinh(u)
            [u,du]=get_props(u);
            val=['sinh(',u,')'];
            der=rise_algebra_cat('/',du,['cosh(',u,')']);
            obj=sad_forward(val,der);
        end
        function obj=tan(u)
            [u,du]=get_props(u);
            val=['tan(',u,')'];
            der=rise_algebra_cat('^',['cos(',u,')'],'2');
            der=rise_algebra_cat('/',du,der);
            obj=sad_forward(val,der);
        end
        function obj=atan(u)
            [u,du]=get_props(u);
            val=['atan(',u,')'];
            der=rise_algebra_cat('^',u,'2');
            der=rise_algebra_cat('+','1',der);
            der=rise_algebra_cat('/',du,['sqrt(',der,')']);
            obj=sad_forward(val,der);
        end
        function obj=tanh(u)
            [u,du]=get_props(u);
            val=['tanh(',u,')'];
            der=rise_algebra_cat('^',['cosh(',u,')'],'2');
            der=rise_algebra_cat('/',du,der);
            obj=sad_forward(val,der);
        end
        function obj=min(u,v)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            val=['min(',u,',',v,')'];
            
            der=rise_algebra_cat('+',...
                rise_algebra_cat('*',uLv,du),...
                rise_algebra_cat('*','(1+uminus(',uLv,'))',dv));
            obj = sad_forward(val,der);
        end
        function obj=max(u,v)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            val=['max(',u,',',v,')'];
            der=rise_algebra_cat('+',rise_algebra_cat('*',uLv,dv),rise_algebra_cat('*','(1+uminus(',uLv,'))',du));
            obj = sad_forward(val,der);
        end
        function obj=sum(u,v)
            if nargin==1
                [u0,du0]=get_props(u(1));
                val=u0;
                der=du0;
                for ii=2:numel(u)
                    [u0,du0]=get_props(u(ii));
                    val=rise_algebra_cat(val,u0,'plus');
                    der=rise_algebra_cat(der,du0,'plus');
                end
                obj=sad_forward(val,der);
            else
                obj=plus(u,v);
            end
        end
        function obj=normpdf(u,mu,sig)
            if nargin<3
                sig='1';
                if nargin<2
                    mu='0';
                end
            end
            [u,du]=get_props(u);
            mu=get_props(mu);
            sig=get_props(sig);
            val=['normpdf(',u,',',mu,',',sig,')'];
            der0=rise_algebra_cat('-',u,mu);
            der1=rise_algebra_cat('^',sig,'2');
            der=rise_algebra_cat('/',['uminus(',der0,')'],der1);
            der=rise_algebra_cat('*',der,du);
            der=rise_algebra_cat('*',der,val);
            obj = sad_forward(val,der);
        end
        function obj=normcdf(u,mu,sig)
            if nargin<3
                sig='1';
                if nargin<2
                    mu='0';
                end
            end
            [u,du]=get_props(u);
            mu=get_props(mu);
            sig=get_props(sig);
            val=['normcdf(',u,',',mu,',',sig,')'];
            der=rise_algebra_cat('*',du,['normpdf(',u,',',mu,',',sig,')']);
            obj = sad_forward(val,der);
        end
        function obj=abs(u)
            [u,du]=get_props(u);
            val=['abs(',u,')'];
            der=rise_algebra_cat('*',['sign(',u,')'],du);
            obj=sad_forward(val,der);
        end
        function obj=isreal(u)
            u=get_props(u);
            val=['real(',u,')'];
            der='0';
            obj=sad_forward(val,der);
        end
        function obj=sqrt(u)
            [u,du]=get_props(u);
            xx=['sqrt(',u,')'];
            dxx=rise_algebra_cat('*','0.5',rise_algebra_cat('/',du,xx));
            obj=sad_forward(xx,dxx);
        end
        function obj=norm(u)
            obj = sqrt(sum(u.^2));
        end
        
        function d = char(u)
            d = u.dx;
        end
    end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
    end
end

function [u,du]=get_props(x)
du='0';
switch class(x)
    case 'sad_forward'
        u=x.x;
        du=x.dx;
    case 'char'
        u=x;
    case 'double'
		u=sprintf('%0.10g',x); % <-- u=num2str(x,10);
    otherwise
        error([mfilename,':: unsupported class ',class(x)])
end
end

