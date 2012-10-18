classdef automatic
    % inspired from autodiff (Matlab central)
    properties
        x
        dx
    end
    methods
        function obj=automatic(a,b)
            % AUTOMATIC automatic class constructor
            %   obj = automatic(a) creates a automatic object with value =a
            %   and derivative =1
            %   obj = automatic(a,b) sets the derivative to b
            if nargin
                na=numel(a);
                if nargin<2
                    b=eye(na);
                end
                obj.x=a; obj.dx=b;
            end
        end
        function y=subsasgn(y,S,x)
            if  isempty(S.subs{1}) 
                return 
            end
            if length(S.subs)>1 
                error('No matrix assignments!'); 
            end
            
            if isa(x, 'automatic')
                if isa(y, 'automatic')
                    y.x(S.subs{:},1) = x.x;
                    y.dx(S.subs{1},:) = x.dx;
                else
                    nd = size(x.dx,2);
                    ny = length(y);
                    testDeriv = zeros(ny,nd);
                    y = automatic(y,testDeriv);
                    y.x(S.subs{:},1) = x.x;
                    y.dx(S.subs{1},:) = x.dx;
                end
            else
                y.x(S.subs{:},1) = x;
                y.dx(S.subs{1},:) = 0;
            end
        end
        function x = subsref(x, S)
            x.x = x.x(S.subs{:});
            x.dx = x.dx(S.subs{1},:);
        end
        function prop=get(x,prop)
            if strcmp(prop,'dx')
                 prop=vertcat(x.(prop));
            else
            prop=x.(prop);
            end
        end
        function x=set(x,prop,val)
            x.(prop)=val;
        end
        function h = normpdf(x,mu,sig)
            % automatic/NORMPDF overloads normpdf with a automatic object argument
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            % simply write the problem in terms of elementary function and
            % let the procedure do the rest
            h=1./(sig.*sqrt(2.*pi)).*exp(-.5.*((x-mu)./sig).^2);
        end
        function h = normcdf(x,mu,sig)
            % automatic/NORMCDF overloads normcdf with a automatic object argument
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            hval=normcdf(x.x,mu,sig);
            der=1./(sig.*sqrt(2.*pi)).*exp(-.5.*((x.x-mu)./sig).^2);
            h = automatic(hval,der);
        end
        function c = mldivide(u,v)
            % automatic/MLDIVIDE overloads mldivide with a automatic object argument
            switch [class(u),class(v)]
                case 'doubleautomatic'
                    c = automatic(u.\v.x, u.\v.dx);
                case {'automaticdouble','automaticautomatic'}
                    if length(u)==1
                        c = times(v,1./u);
                    else
                        error('Cannot do u.\v when u is an automatic object with length>1')
                    end
                otherwise
                    error(['.\ not defined for ',class(u),' and ',class(v)]);
            end
        end
        function c = minus(u,v)
            % automatic/MINUS overloads minus with a automatic object argument
            c = plus(u,-v);
        end
        function c = min(u,v)
            % automatic/MIN overloads min with a automatic object argument
            % min(u) selects the value of u which is u minimum.
            % Both min(u) and min(u,v) will work.
            
            if nargin==1
                [y,i] = min(u.x);
                c = automatic(y,u(i).dx);
            else
                switch [class(u),class(v)]
                    case 'doubleautomatic'
                        ok = v.x<u;
                        c = automatic(min(u,v.x), ok.*v.dx);
                    case 'automaticdouble'
                        ok = u.x<v;
                        c = automatic(min(u.x,v), ok.*u.dx);
                    case 'automaticautomatic'
                        ok = u.x<v.x;
                        c = automatic(min(u.x,v.x), ok.*u.dx+(1-ok).*v.dx);
                    otherwise
                        error(['Can''t compute min of ',class(u),' and ',class(v)]);
                end
            end
        end
        function c = max(u,v)
            % automatic/MAX overloads max with a automatic object argument
            % max(u) selects the value of u which is u maximum.
            % Both max(u) and max(u,v) will work.
            
            if nargin==1
                [y,i] = max(u.x);
                c = automatic(y,u(i).dx);
            else
                switch [class(u),class(v)]
                    case 'doubleautomatic'
                        ok = v.x>u;
                        c = automatic(min(u,v.x), ok.*v.dx);
                    case 'automaticdouble'
                        ok = u.x>v;
                        c = automatic(min(u.x,v), ok.*u.dx);
                    case 'automaticautomatic'
                        ok = u.x>v.x;
                        c = automatic(min(u.x,v.x), ok.*u.dx+(1-ok).*v.dx);
                    otherwise
                        error(['Can''t compute max of ',class(u),' and ',class(v)]);
                end
            end
        end
        function c = log10(u)
            % automatic/LOG10 overloads log10 with a automatic object argument
            c = automatic(log(u.x)./log(10), (1./u.x./log(10).*u.dx));
        end
        function c = log(u)
            % automatic/LOG overloads log with a automatic object argument
            c = automatic(log(u.x),(1./u.x.*u.dx));
        end
        function c = isreal(u)
            % automatic/ISREAL overloads isreal with a automatic object argument
            c = isreal(u.x);
        end
        function c = exp(u)
            % automatic/EXP overloads exp with a automatic object argument
            c = automatic(exp(u.x),(exp(u.x).*u.dx));
        end
        function d = double(u)
            % automatic/DOUBLE overloads double with a automatic object
            % argument. Returns the value and the derivatives concatenated
            d = [vertcat(u.x),vertcat(u.dx)];
        end
        function da = diff(ad,n)
            % automatic/DIFF overloads diff with a automatic object argument
            if nargin<2
                n=1;
            end
            da = automatic(diff(ad.x,n),diff(ad.dx,n));
        end
        function c = cosh(u)
            % automatic/COSH overloads cosh with a automatic object argument
            c = automatic(cosh(u.x),(sinh(u.x).*u.dx));
        end
        function c = cos(u)
            % automatic/COS overloads cos with a automatic object argument
            c = automatic(cos(u.x),(-sin(u.x).*u.dx));
        end
        function c = conj(u)
            % automatic/CONJ overloads conj with a automatic object argument
            c=automatic(conj(u.x),conj(u.dx));
        end
        function c = atan(u)
            % automatic/ATAN overloads atan with a automatic object argument
            c = automatic(atan(u.x),(1./sqrt(1+u.x.^2).*u.dx));
        end
        function c = asin(u)
            % automatic/ASIN overloads asin with a automatic object argument
            c = automatic(asin(u.x),(1./sqrt(1-u.x.^2).*u.dx));
        end
        function c = acos(u)
            % automatic/ACOS overloads acos with a automatic object argument
            c = automatic(acos(u.x),(-1./sqrt(1-u.x.^2).*u.dx));
        end
        function c = abs(ad)
            % automatic/ABS overloads abs with a automatic object argument
            if all(isreal(ad.x))&& all(isreal(ad.dx(:)))
                c = automatic(abs(ad.x),(sign(ad.x).*ad.dx));
            else
                c = sqrt(ad.*conj(ad));
            end
        end
        function c = uminus(ad)
            % automatic/UMINUS overloads uminus with a automatic object argument
            c = automatic( -ad.x, -ad.dx);
        end
        function c = times(u,v)
            % automatic/TIMES overloads times with a automatic object argument
            switch [class(u),class(v)]
                case 'automaticdouble'
                    c = automatic(v.*u.x, v.*u.dx);
                case 'doubleautomatic'
                    c = automatic(u.*v.x, u.*v.dx);
                case 'automaticautomatic'
                    c = automatic( u.x.*v.x, v.x.*u.dx+u.x.*v.dx);
                otherwise
                    error(['Can''t multiply (.*) ',class(u),' and ',class(v)]);
            end
        end
        function c = tanh(u)
            % automatic/TANH overloads tanh with a automatic object argument
            c = automatic(tanh(u.x), 1./cosh(u.x).^2.*u.dx);
        end
        function c = tan(u)
            % automatic/TAN overloads tan with a automatic object argument
            c = automatic(tan(u.x), 1./cos(u.x).^2.*u.dx);
        end
        function c = sum(u,v)
            % automatic/SUM overloads sum with a automatic object argument
            if nargin==1
                c = automatic(sum([u.x]), sum([u.dx]));
            else
                c=plus(u,v);
            end
        end
        function c = sqrt(u)
            % automatic/SQRT overloads sqrt with a automatic object argument
            c = u.^(0.5);
        end
        function c = sinh(u)
            % automatic/SINH overloads sinh with a automatic object argument
            c = automatic( sinh(u.x), (cosh(u.x).*u.dx));
        end
        function c = real(u)
            % automatic/REAL overloads real with a automatic object argument
            c=automatic(real(u.x),real(u.dx));
        end
        function c = rdivide(u,v)
            % automatic/RDIVIDE overloads rdivide with a automatic object argument
            switch [class(u),class(v)]
                case 'automaticdouble'
                    c = times(u,1./v);
                    
                case {'doubleautomatic','automaticautomatic'}
                    c = times(u,reciprocal(v));
                    
                otherwise
                    error(['Can''t divide (./) ',class(u),' and ',class(v)]);
            end
        end
        function c = power(u,v)
            % automatic/POWER overloads power with a automatic object argument
            switch [class(u),class(v)]
                case 'automaticdouble'
                    c = automatic( u.x.^v, (v.*(u.x.^(v-1)).*u.dx));
                case 'doubleautomatic'
                    c = automatic( u.^v.x, u.^v.x.*log(u).*v.dx);
                case 'automaticautomatic'
                    c = exp(log(u).*v);
                otherwise
                    error(['Can''t do ',class(u),'.^',class(v)]);
            end
        end
        function c = plus(u,v)
            % automatic/PLUS overloads  with a automatic object argument
            switch [class(u),class(v)]
                case 'automaticdouble'
                    c = automatic(u.x+v, u.dx);
                case 'doubleautomatic'
                    c = automatic(v.x+u, v.dx);
                case 'automaticautomatic'
                    c = automatic(u.x+v.x,u.dx+v.dx);
                otherwise
                    error(['Can''t add ',class(u),' and ',class(v)]);
            end
        end
        function n = norm(u)
            % automatic/NORM overloads norm with a automatic object argument
            n = sqrt(sum(u.^2));
        end
        function c = mtimes(u,v)
            % automatic/MTIMES overloads mtimes with a automatic object argument
            if all(size(u)==1) || all(size(v)==1)
                c = times(u,v);
                return
            end
            switch [class(u),class(v)]
                case 'automaticdouble'
                    c = automatic(u.x.*v,v.*u.dx);
                case 'doubleautomatic'
                    c = automatic(u.*v.x,u.*v.dx);
                case 'automaticautomatic'
                    c = automatic(u.x.*v.x, u.x.*v.dx+v.x.*u.dx);
                otherwise
                    error(['Can''t matrix multiply ',class(u),' and ',class(v)]);
                    
            end
        end
        function c = mrdivide(u,v)
            % automatic/MRDIVIDE overloads mrdivide with a automatic object argument
            if numel(v)==1
                c = times(u,1./v);
            else
                error('Cannot use ./ for automatic objects when denominator is nonscalar');
            end
        end
        function c = mpower(u,v)
            % automatic/MPOWER overloads mpower with a automatic object argument
            c=power(u,v);
        end
    end
    methods(Access=private)
        function ad1 = reciprocal(ad)
            % RECIPROCAL computes the reciprocal of an automatic object,
            % used for division.
            ad1 = automatic(1./ad.x,-1./ad.x.^2.*ad.dx);
        end
    end
end

