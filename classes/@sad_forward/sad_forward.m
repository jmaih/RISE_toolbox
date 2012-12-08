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
                if ~ischar(y)
                    error('first input must be char')
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
            xx='0';
            dxx='0';
            if ~strcmp(u,'0')
                xx=strcat('-',parenthesize(u,'-'));
            end
            if ~strcmp(du,'0')
                dxx=strcat('-',parenthesize(du,'-'));
            end
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
            dxx1=rise_algebra_cat('*',rise_algebra_cat('*',dv,strcat('log(',u,')')),xx);
            dxx2=rise_algebra_cat('*',rise_algebra_cat('*',v,du),...
                rise_algebra_cat('^',u,rise_algebra_cat('-',v,'1')));
            dxx=rise_algebra_cat('+',dxx1,dxx2);
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
            xx=strcat('exp(',u,')');
            dxx=rise_algebra_cat('*',du,xx);
            obj=sad_forward(xx,dxx);
        end
        function obj=log(u)
            [u,du]=get_props(u);
            xx=strcat('log(',u,')');
            dxx=rise_algebra_cat('/',du,u);
            obj=sad_forward(xx,dxx);
        end
        function obj=log10(u)
            obj=log(u)/log(10);
        end
        function obj=cos(u)
            [u,du]=get_props(u);
            xx=strcat('cos(',u,')');
            dxx=rise_algebra_cat('*',...
                rise_algebra_cat('-','0',du),strcat('sin(',u,')'));
            obj=sad_forward(xx,dxx);
        end
        function obj=acos(u)
            [u,du]=get_props(u);
            val=strcat('acos(',u,')');
            der=rise_algebra_cat('^',u,'2');
            der=rise_algebra_cat('-','1',der);
            der=rise_algebra_cat('/',rise_algebra_cat('*','-1',du),...
                strcat('sqrt(',der,')'));
            obj=sad_forward(val,der);
        end
        function obj=cosh(u)
            [u,du]=get_props(u);
            val=strcat('cosh(',u,')');
            der=rise_algebra_cat('*',du,strcat('sinh(',u,')'));
            obj=sad_forward(val,der);
        end
        function obj=sin(u)
            [u,du]=get_props(u);
            xx=strcat('sin(',u,')');
            dxx=rise_algebra_cat('*',du,strcat('cos(',u,')'));
            obj=sad_forward(xx,dxx);
        end
        function obj=asin(u)
            [u,du]=get_props(u);
            val=strcat('asin(',u,')');
            der=rise_algebra_cat('^',u,'2');
            der=rise_algebra_cat('*','-1',der);
            der=rise_algebra_cat('/',du,strcat('sqrt(',der,')'));
            obj=sad_forward(val,der);
        end
        function obj=sinh(u)
            [u,du]=get_props(u);
            val=strcat('sinh(',u,')');
            der=rise_algebra_cat('/',du,strcat('cosh(',u,')'));
            obj=sad_forward(val,der);
        end
        function obj=tan(u)
            [u,du]=get_props(u);
            val=strcat('tan(',u,')');
            der=rise_algebra_cat('^',strcat('cos(',u,')'),'2');
            der=rise_algebra_cat('/',du,der);
            obj=sad_forward(val,der);
        end
        function obj=atan(u)
            [u,du]=get_props(u);
            val=strcat('atan(',u,')');
            der=rise_algebra_cat('^',u,'2');
            der=rise_algebra_cat('+','1',der);
            der=rise_algebra_cat('/',du,strcat('sqrt(',der,')'));
            obj=sad_forward(val,der);
        end
        function obj=tanh(u)
            [u,du]=get_props(u);
            val=strcat('tanh(',u,')');
            der=rise_algebra_cat('^',strcat('cosh(',u,')'),'2');
            der=rise_algebra_cat('/',du,der);
            obj=sad_forward(val,der);
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
            val=strcat('normpdf(',u,',',mu,',',sig,')');
            der0=rise_algebra_cat('-',u,mu);
            der1=rise_algebra_cat('^',sig,'2');
            der=rise_algebra_cat('/',strcat('uminus(',der0,')'),der1);
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
            val=strcat('normcdf(',u,',',mu,',',sig,')');
            der=rise_algebra_cat('*',du,strcat('normpdf(',u,',',mu,',',sig,')'));
            obj = sad_forward(val,der);
        end
        function obj=abs(u)
            [u,du]=get_props(u);
            val=strcat('abs(',u,')');
            der=rise_algebra_cat('*',strcat('sign(',u,')'),du);
            obj=sad_forward(val,der);
        end
        function obj=sqrt(u)
            [u,du]=get_props(u);
            xx=strcat('sqrt(',u,')');
            dxx=rise_algebra_cat('*','0.5',rise_algebra_cat('/',du,xx));
            obj=sad_forward(xx,dxx);
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

