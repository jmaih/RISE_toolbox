classdef rise_sad
    % Symbolic automatic differentiation
    % inspired from autodiff (Matlab central)
    properties
        x
        dx
    end
    methods
        function obj=rise_sad(a,b)
            %  rise_sad class constructor
            %   obj = rise_sad(a) creates a rise_sad object with value =a
            %   and derivative =1
            %   obj = rise_sad(a,b) sets the derivative to b
            if nargin
                if isa(a,'rise_sad')
                    obj=a;
                else
                    if ischar(a),a=cellstr(a); end
                    na=numel(a);
                    if nargin<2
                        b=repmat({'1'},na,1);
                    end
                    if ischar(b),b=cellstr(b); end
                    obj=rise_sad.empty(0,1);
                    for ii=1:na
                        %                        obj(ii).x={Optimize(a{ii})};  % Optimize();
                        %                        obj(ii).dx={Optimize(b{ii})}; % Optimize()
                        obj(ii).x=a(ii);  % Optimize();
                        obj(ii).dx=b(ii); % Optimize()
                    end
                end
            end
        end
        %---------- functions in two arguments ------------%
        function c = plus(u,v)
            % rise_sad/PLUS overloads  with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uPv=myplus(u,v);
            duPdv=myplus(du,dv);
            c = rise_sad(uPv,duPdv);
        end
        function c = minus(u,v)
            % rise_sad/MINUS overloads minus with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uMv=myminus(u,v);
            duMdv=myminus(du,dv);
            c = rise_sad(uMv,duMdv);
        end
        function c = mldivide(u,v)
            % rise_sad/MLDIVIDE overloads mldivide with a rise_sad object argument
            c = mrdivide(v,u);
        end
        function c = mrdivide(u,v)
            % rise_sad/MRDIVIDE overloads mrdivide with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            hval=mydivide(u,v);
            if strcmp(dv,'0')
                der=mydivide(du,v);
            else
                der=mytimes(du,v);
                der=myminus(der,mytimes(dv,u));
                der=mydivide(der,mypower(v,'2'));
            end
            c = rise_sad(hval,der);
            %             c = rise_sad(['(',u,')/(',v,')'],['((',du,')*(',v,')-(',dv,')*(',u,'))/(',v,')^2']);
        end
        function c = rdivide(u,v) %
            % rise_sad/RDIVIDE overloads rdivide with a rise_sad object argument
            c = mrdivide(u,v);
        end
        function c = ldivide(u,v) %
            % rise_sad/LDIVIDE overloads ldivide with a rise_sad object argument
            c = rdivide(v,u);
        end
        function c = times(u,v)
            % rise_sad/TIMES overloads times with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            der=[parenthesize(du),'.*',parenthesize(v),'+',parenthesize(dv),'.*',parenthesize(u)];
            c = rise_sad([parenthesize(u),'.*',parenthesize(v)],der);
        end
        function c = mtimes(u,v)
            % rise_sad/MTIMES overloads mtimes with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            hval=mytimes(u,v);
            der=mytimes(du,v);
            der=myplus(der,mytimes(dv,u));
            c = rise_sad(hval,der);
%             c = rise_sad(['(',u,')*(',v,')'],['(',du,')*(',v,')+(',dv,')*(',u,')']);
        end
        function c = power(u,v)
            % rise_sad/POWER overloads power with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uPv=[parenthesize(u),'.^',parenthesize(v)];
            der=['(',parenthesize(dv),'.*log(',u,')+',parenthesize(du),'./',parenthesize(u),'.*',parenthesize(v),').*',uPv];
            c = rise_sad(uPv,der);
        end
        function c = mpower(u,v)
            % rise_sad/MPOWER overloads mpower with a rise_sad object argument
            [u,du]=get_props(u);[v,dv]=get_props(v);
            hval=mypower(u,v);
            if strcmp(dv,'0')
                der=mytimes(mytimes(v,du),mypower(u,myminus(v,'1')));
            else
                der1=mytimes(dv,['log(',u,')']);
                der2=mytimes(mydivide(du,u),v);
                der=mytimes(myplus(der1,der2),hval);
            end
            c = rise_sad(hval,der);
            %             uPv=['(',u,')^(',v,')'];
%             c = rise_sad(uPv,['((',dv,')*log(',u,')+(',du,')/(',u,')*(',v,'))*',uPv]);
        end
        function c = min(u,v)
            % rise_sad/MIN overloads min with a rise_sad object argument
            % but will work with 2 arguments.
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            c=rise_sad(['min(',u,',',v,')'],[uLv,'*',parenthesize(du),'+(1-',uLv,')*',parenthesize(dv)]);
        end
        function c = max(u,v)
            % rise_sad/MAX overloads max with a rise_sad object argument
            % max(u) selects the value of u which is u maximum.
            % Both max(u) and max(u,v) will work.
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            c=rise_sad(['max(',u,',',v,')'],[uLv,'*',parenthesize(dv),'+(1-',uLv,')*',parenthesize(du)]);
        end
        %---------- functions in one argument ------------%
        function c = exp(u)
            % rise_sad/EXP overloads exp with a rise_sad object argument
            [u,du]=get_props(u);
            hval=['exp(',u,')'];
            der=mytimes(du,hval);
            c = rise_sad(hval,der);
        end
        function c = log10(u)
            % rise_sad/LOG10 overloads log10 with a rise_sad object argument
            c=rdivide(log(u),'log(10)');
        end
        function c = log(u)
            % rise_sad/LOG overloads log with a rise_sad object argument
            [u,du]=get_props(u);
            c = rise_sad(['log(',u,')'],mydivide(du,u));
%             c = rise_sad(['log(',u,')'], ['(',du,')/(',u,')']);
        end
        function c = isreal(u)
            % rise_sad/ISREAL overloads isreal with a rise_sad object argument
            u=get_props(u);
            c = ['isreal(',u,')'];
        end
        function d = char(u)
            % rise_sad/DOUBLE overloads double with a rise_sad object
            % argument. Returns the value and the derivatives concatenated
            d = u.dx{1};
        end
        function c = cosh(u)
            % rise_sad/COSH overloads cosh with a rise_sad object argument
            [u,du]=get_props(u);
            der=mytimes(du,['sinh(',u,')']);
            c = rise_sad(['cosh(',u,')'],der);
        end
        function c = cos(u)
            % rise_sad/COS overloads cos with a rise_sad object argument
            [u,du]=get_props(u);
            der=mytimes(['-sin(',u,')'],du);
            c = rise_sad(['cos(',u,')'],der);
        end
        function c = conj(u)
            % rise_sad/CONJ overloads conj with a rise_sad object argument
            [u,du]=get_props(u);
            c=rise_sad(['conj(',u,')'],['conj(',du,')']);
        end
        function c = atan(u)
            % rise_sad/ATAN overloads atan with a rise_sad object argument
            [u,du]=get_props(u);
            der=mypower(u,'2');
            der=myplus('1',der);
            der=mydivide(du,['sqrt(',der,')']);
            c = rise_sad(['atan(',u,')'],der);
%             c = rise_sad(['atan(',u,')'],['(',du,')/sqrt(1+(',u,')^2)']);
        end
        function c = asin(u)
            % rise_sad/ASIN overloads asin with a rise_sad object argument
            [u,du]=get_props(u);
            der=mypower(u,'2');
            der=myminus('1',der);
            der=mydivide(du,['sqrt(',der,')']);
            c = rise_sad(['asin(',u,')'],der);
%             c = rise_sad(['asin(',u,')'],['(',du,')/sqrt(1-(',u,')^2)']);
        end
        function c = acos(u)
            % rise_sad/ACOS overloads acos with a rise_sad object argument
            [u,du]=get_props(u);
            der=mypower(u,'2');
            der=myminus('1',der);
            der=mydivide(['-',parenthesize(du)],['sqrt(',der,')']);
            c = rise_sad(['acos(',u,')'],der);
%             c = rise_sad(['acos(',u,')'],['-(',du,')/sqrt(1-(',u,')^2)']);
        end
        function c = abs(u)
            % rise_sad/ABS overloads abs with a rise_sad object argument.
            % but works only for real variables
            [u,du]=get_props(u);
            c = rise_sad(['abs(',u,')'],mytimes(['sign(',u,')'],du));
        end
        function c = uminus(u)
            % rise_sad/UMINUS overloads uminus with a rise_sad object argument
            c = rise_sad(['-',parenthesize(u.x{1})], ['-',parenthesize(u.dx{1})]);
        end
        function c = uplus(u)
            % rise_sad/UMINUS overloads uminus with a rise_sad object argument
            c = rise_sad(u.x{1},u.dx{1});
        end
        function c = tanh(u)
            % rise_sad/TANH overloads tanh with a rise_sad object argument
            [u,du]=get_props(u);
            der=mypower(['cosh(',u,')'],'2');
            der=mydivide(du,der);
            c = rise_sad(['tanh(',u,')'],der);
        end
        function c = tan(u)
            % rise_sad/TAN overloads tan with a rise_sad object argument
            [u,du]=get_props(u);
            der=mypower(['cos(',u,')'],'2');
            der=mydivide(du,der);
            c = rise_sad(['tan(',u,')'],der);
        end
        function c = sum(u,v)
            % rise_sad/SUM overloads sum with a rise_sad object argument
            if nargin==1
                [u0,du0]=get_props(u(1));
                xx=u0;
                dxx=du0;
                for ii=2:numel(u)
                    [u0,du0]=get_props(u(ii));
                    xx=myplus(xx,u0); 
                    dxx=myplus(dxx,du0);
                end
                c = rise_sad(xx,dxx);
            else
                c=plus(u,v);
            end
        end
        function c = sqrt(u)
            % rise_sad/SQRT overloads sqrt with a rise_sad object argument
            c = u^(0.5);
        end
        function c = sinh(u)
            % rise_sad/SINH overloads sinh with a rise_sad object argument
            [u,du]=get_props(u);
            der=mydivide(du,['cosh(',u,')']);
            c = rise_sad(['sinh(',u,')'],der);
        end
        function c = real(u)
            % rise_sad/REAL overloads real with a rise_sad object argument
            [u,du]=get_props(u);
            c=rise_sad(['real(',u,')'],['real(',du,')']);
        end
        function n = norm(u)
            % rise_sad/NORM overloads norm with a rise_sad object argument
            n = sqrt(sum(u.^2));
        end
        %---------- functions in one argument and varargin ------------%
        function h = normpdf(u,mu,sig)
            % rise_sad/NORMPDF overloads normpdf with a rise_sad object argument
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            [u,du]=get_props(u);
            [mu,sig]=get_char_form(mu,sig);
            % simply write the problem in terms of elementary function and
            % let the procedure do the rest
            hval=['normpdf(',u,',',mu,',',sig,')'];
            der0=myminus(u,mu);
            der1=mypower(sig,'2');
            der=mydivide(['-(',der0,')'],der1);
            der=mytimes(der,du);
            der=mytimes(der,hval);
% %             der=['-(',u,'-(',mu,'))/(',sig,')^2*(',du,')*',hval];
            h = rise_sad(hval,der);
        end
        function h = normcdf(u,mu,sig)
            % rise_sad/NORMCDF overloads normcdf with a rise_sad object argument
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            [u,du]=get_props(u);
            [mu,sig]=get_char_form(mu,sig);
            hval=['normcdf(',u,',',mu,',',sig,')']; %<--- hval=['0.5*(1+erf((',x.x{1},'-',mu,')/(',sig,'*sqrt(2*pi))))']; % <--- hval=['normcdf(',x.x{1},',',mu,',',sig,')'];
            der=mytimes(du,['normpdf(',u,',',mu,',',sig,')']);
            h = rise_sad(hval,der);
        end
    end
    methods(Static)
        varargout=optimize(varargin)
        varargout=diff(varargin)
        varargout=hessian(varargin)
        varargout=jacobian(varargin)
    end
end

function varargout=get_char_form(varargin)
n=length(varargin);
for ivar=1:n
    vv=varargin{ivar};
    switch class(vv)
        case 'char'
            varargout{ivar}=vv;
        case 'double'
            varargout{ivar}=num2str(vv,10);
        case 'rise_sad'
            varargout{ivar}=vv.x{1};
        otherwise
    end
end
end

function [u,du]=get_props(x)
du='0';
switch class(x)
    case 'rise_sad'
        u=x.x{1};
        du=x.dx{1};
    case 'char'
        u=x;
    case 'double'
        u=num2str(x,10);
    otherwise
        error([mfilename,':: unsupported class ',class(x)])
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
        c=tryevaluate(['-',parenthesize(b)]);
    end
else
    if strcmp(b,'0')
        c=a;
    else
        c=tryevaluate([a,'-',parenthesize(b)]);
    end
end
end

function c=mytimes(a,b)
if strcmp(a,'0')||strcmp(b,'0')
    c='0';
else
    c=tryevaluate([parenthesize(a),'*',parenthesize(b)]);
end
end

function c=mypower(a,b)
if strcmp(b,'0')
    c='1';
else
    c=tryevaluate([parenthesize(a),'^',parenthesize(b)]);
end
end

function c=mydivide(a,b)
if strcmp(a,'0')
    c='0';
else
    c=tryevaluate([parenthesize(a),'/',parenthesize(b)]);
end
end

function x=parenthesize(x)
flag=false;
for ii=1:length(x)
    if any(x(ii)=='+-*/^')
        flag=true;
        break
    end
end
if flag
    x=['(',x,')'];
end
end


