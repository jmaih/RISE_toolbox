classdef dyn_series
    properties
        val
        coef
    end
    methods
        function obj = dyn_series(a,c,order)
            %DYN_SERIES dyn_series class constructor
            %   Modifies Richard Neidinger's series
            %   obj = dyn_series(a,c,order) creates a dyn_series object for a Talyor polynomial with
            %   constant term a, following coefficient(s) in c, and zeros thru degree order
            %   put into fields:
            %       .val for the constant term value of the function
            %       .coef for the vector of coefficients from linear term to highest term
            %   h = dyn_series(a,c) will assume that c contains all desired coefficients
            %   obj = dyn_series(a,1,order) creates a dyn_series for a variable obj at value a.
            if nargin > 0
                if strcmp(class(a),'dyn_series') % isa(a,'dyn_series')
                    obj = a;
                else
                    obj.val = a;
                    if nargin>1
                        obj.coef = c;
                        if nargin >2
                            if numel(c)<order
                                obj.coef = [obj.coef zeros(1,order-length(c))];
                            elseif numel(c)>order
                                error([mfilename,':: the second input cannot have more elements than the order of approximation(',int2str(order),')'])
                            end
                        end
                    end
                end
            end
        end
        %=============================
        function obj = assigncoef1(obj,w)
            %ASSIGNCOEF1 assigns coef(1) of dyn_series object obj to value in w.
            %   obj may be a vector of dyn_series objects and w a corresponding vector of values.
            for i = 1:numel(obj)
                obj(i).coef(1) = w(i);
            end
        end
        %=============================
        function vec = coefs(obj)
            %COEFS takes dyn_series object and returns vector of coefficients from linear to highest term.
            vec = obj.coef;
        end
        %=============================
        function coefs = double(f)
            %dyn_series/DOUBLE Convert dyn_series object to vector of doubles
            %   containing all the dyn_series coefficient values.
            coefs = [vertcat(f.val) vertcat(f.coef)];
        end
        %=============================
        function h = minus(u,v)
            %dyn_series/MINUS overloads subtraction with at least one dyn_series object argument
            % assumes .coef lists are of same length
            type =isa(u,'dyn_series')+2*isa(v,'dyn_series');
            switch type
                case 1
                    h=dyn_series(u.val-v,u.coef);
                case 2
                    h=dyn_series(u-v.val,-v.coef);
                case 3
                    h=dyn_series(u.val-v.val,u.coef-v.coef);
            end
        end
        %=============================
        function h = uminus(u)
            %dyn_series/UMINUS overloads negation with a dyn_series object argument
            h = dyn_series(-u.val,-u.coef);
        end
        %=============================
        function h = plus(u,v)
            %dyn_series/PLUS overloads addition with at least one dyn_series object argument
            % assumes .coef lists are of same length
            type =isa(u,'dyn_series')+2*isa(v,'dyn_series');
            switch type
                case 1
                    h=dyn_series(u.val+v,u.coef);
                case 2
                    h=dyn_series(u+v.val,v.coef);
                case 3
                    h=dyn_series(u.val+v.val,u.coef+v.coef);
            end
        end
        %=============================
        function h = sqrt(u)
            %dyn_series/SQRT overloads square root of a dyn_series object argument
            order = length(u.coef);
            hval=sqrt(u.val);
            h02i=1/(2*hval);
            hcoef = zeros(1,order);
            for ii=1:order
                hcoef(ii)=h02i*(u.coef(ii)-dot(hcoef(1:ii-1),hcoef(ii-1:-1:1)));
            end
            h=dyn_series(hval,hcoef);
        end
        %=============================
        function h = exp(u)
            %dyn_series/EXP overloads exp with a dyn_series object argument
            order = length(u.coef);
            up = (1:order).*u.coef;
            hvec = [exp(u.val),zeros(1,order)];
            for k = 1:order
                hvec(k+1) = dot(up(1:k),hvec(k:-1:1))/k;
            end
            h = dyn_series(hvec(1), hvec(2:order+1));
        end
        %=============================
        function h = log(u)
            %dyn_series/LOG overloads natural logarithm of a dyn_series object argument
            order = length(u.coef);
            hp = [u.coef(1)/u.val,zeros(1,max(0,order-1))];
            for k = 2:order
                hp(k) = (u.coef(k)-dot((1:k-1).*hp(1:k-1),u.coef(k-1:-1:1))/k)/u.val;
            end
            h = dyn_series(log(u.val),hp);
        end
        %=============================
        function h = normpdf(x,mu,sig)
            %dyn_series/NORMPDF overloads normpdf with a dyn_series object argument
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            % put the problem under the form h'=u'*v. In this case, v turns
            % out to be h. So we have h'=u'*h
            uu=-.5*((x-mu)/sig)^2;
            hval=normpdf(x.val,mu,sig);
            order = length(x.coef);
            uup = (1:order).*uu.coef;
            hvec = [hval,zeros(1,order)];
            for k = 1:order
                hvec(k+1) = dot(uup(1:k),hvec(k:-1:1))/k;
            end
            h = dyn_series(hvec(1), hvec(2:order+1));
        end
        %=============================
        function h = normcdf(x,mu,sig)
            %dyn_series/NORMCDF overloads normcdf with a dyn_series object argument
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            hval=normcdf(x.val,mu,sig);
            order = length(x.coef);
            % put the problem under the form h'=u'*v. In this case, we let
            % v=normpdf and u'=1, implying u=x, which we already have
            v=normpdf(x,mu,sig);
            vcoefs=double(v);
            uup = (1:order).*x.coef;
            hvec = [hval,zeros(1,order)];
            for k = 1:order
                hvec(k+1) = dot(uup(1:k),vcoefs(k:-1:1))/k;
            end
            h = dyn_series(hvec(1), hvec(2:order+1));
        end
        %=============================
        function h = power(u,r)
            %dyn_series/POWER overloads power with at least one dyn_series object argument.
            % Assumes .coef lists are of same length, if both are dyn_series objects.
            h = exp(r*log(u));  % if u is scalar, only exp is dyn_series "convolution"
        end
        %=============================
        function h = mpower(u,r)
            %dyn_series/MPOWER overloads power with at least one dyn_series object argument.
            % Assumes .coef lists are of same length, if both are dyn_series objects.
            if isa(r,'dyn_series')
                h = exp(r*log(u));  % if u is scalar, only exp is dyn_series "convolution"
                % if u and r are dyn_series, exp, log, and * are "convolutions"
            else % assume r is a scalar
                if r == 0
                    h = 1;
                elseif r == 1
                    h = u;
                elseif r == 2
                    h = u*u;
                elseif r == 1/2
                    h = sqrt(u);
                elseif u.val==0 && round(r)==r % make this case more efficient later
                    h = u;
                    for k = 2:r
                        h = h*u;
                    end
                else % about same computation as exp&log but this allows negative base
                    order = length(u.coef);
                    up = (1:order).*u.coef;
                    hvec = [u.val^r,zeros(1,order)];
                    hp = zeros(1,order);
                    for k = 1:order
                        hp(k) = (r*dot(hvec(1:k),up(k:-1:1)) ...
                            - dot(hp(1:k-1),u.coef(k-1:-1:1)))/u.val;
                        hvec(k+1) = (hp(k)/k);
                    end
                    h = dyn_series(hvec(1), hvec(2:order+1));
                end
            end
        end
        %=============================
        function h = mrdivide(u,v)
            %dyn_series/MRDIVIDE overloads division with at least one dyn_series object argument
            % assumes .coef lists are of same length
            type =isa(u,'dyn_series')+2*isa(v,'dyn_series');
            switch type
                case 1
                    h=dyn_series(u.val./v,u.coef./v);
                otherwise
                    order = length(v.coef);
                    if type==2
                        % u is a scalar, make it into a series
                        u=dyn_series(u,0,order);
                    end
                    u_valcoef=double(u);
                    v_valcoef=double(v);
                    order=numel(u_valcoef)-1;
                    h_coefs=nan(1,order+1);
                    h_coefs(1)=u.val/v.val;
                    v0=v_valcoef(1);
                    for ii=1:order
                        % flipping one of the terms
                        h_coefs(ii+1)=(u_valcoef(ii+1)...
                            -dot(h_coefs(1:ii),v_valcoef(ii+1:-1:2))...
                            )/v0;
                    end
                    h=dyn_series(h_coefs(1),h_coefs(1+(1:order)));
            end
        end
        %=============================
        function h = times(u,v)
            %dyn_series/TIMES overloads multiplication with at least one dyn_series object argument
            % assumes .coef lists are of same length
            usiz=size(u);
            vsiz=size(v);
            if isequal(usiz,vsiz)
                h(usiz)=dyn_series;
                for ii=1:prod(usiz)
                    h(ii)=u(ii)*v(ii);
                end
            elseif max(vsiz)>1 && max(usiz)==1
                q=num2cell(vsiz);
                h(q{:})=dyn_series;
                for ii=1:prod(vsiz)
                    h(ii)=u*v(ii);
                end
            elseif max(vsiz)==1 && max(usiz)>1
                q=num2cell(usiz);
                h(q{:})=dyn_series;
                for ii=1:prod(usiz)
                    h(ii)=u(ii)*v;
                end
            else
                error([mfilename,':: elements must have the same size or at least one of them should be a scalar'])
            end
        end
        %=============================
        function h = mtimes(u,v)
            %dyn_series/MTIMES overloads multiplication with at least one dyn_series object argument
            % assumes .coef lists are of same length
            type =isa(u,'dyn_series')+2*isa(v,'dyn_series');
            switch type
                case 1
                    h=dyn_series(u.val*v,u.coef*v);
                case 2
                    h=dyn_series(u*v.val,u*v.coef);
                case 3
                    u_valcoef=double(u);
                    v_valcoef=double(v);
                    order=numel(u_valcoef)-1;
                    h_coefs=nan(1,order);
                    for ii=1:order
                        % flipping one of the terms
                        h_coefs(ii)=dot(u_valcoef(1:ii+1),v_valcoef(ii+1:-1:1));
                    end
                    h=dyn_series(u.val*v.val,h_coefs);
            end
        end
        %=============================
        function g = cos(u)
            %dyn_series/COS overloads cosine of a dyn_series object argument
            [~,g] = sincos(u);
        end
        %=============================
        function h = sin(u)
            %dyn_series/SIN overloads sine with a dyn_series object argument
            h = sincos(u);
        end
        %=============================
        function h = asin(u)
            %dyn_series/ASIN overloads arcsine of a dyn_series object argument
            order = length(u.coef);
            h0 = asin(u.val);
            v0 = cos(h0);
            hp = zeros(1,order);
            vcoeffliplr = zeros(1,order);
            ufliplr = [fliplr(u.coef) u.val];
            for k = 1:order
                hp(k) = (k*u.coef(k) ...
                    - dot(hp(1:k-1),vcoeffliplr(order-(k-2):order))...
                    )/v0;
                vpk = -dot(hp(1:k),ufliplr(order+2-k:order+1));
                vcoeffliplr(order-k+1) = vpk/k;
            end
            h = dyn_series(h0,hp./(1:order));
        end
        %=============================
        function [s, c] = sincos(u)
            %dyn_series/SINCOS overloads exp with a dyn_series object argument
            order = length(u.coef);
            up = (1:order).*u.coef;
            sinvec = sin(u.val);
            cosvec = cos(u.val);
            for k = 1:order
                sinvec(k+1) = dot( cosvec(1:k), up(k:-1:1) ) / k;
                cosvec(k+1) = -dot( sinvec(1:k), up(k:-1:1) ) / k;
            end
            s = dyn_series(sinvec(1), sinvec(2:order+1));
            c = dyn_series(cosvec(1), cosvec(2:order+1));
        end
        %=============================
        function h = tan(u)
            %dyn_series/TAN overloads tangent of a dyn_series object argument
            order = length(u.coef);
            up = (1:order).*u.coef;
            hvec = [tan(u.val),zeros(1,order)];
            hpfliplr = zeros(1,order);
            vAll = [(cos(u.val))^(-2),zeros(1,order)];  % sec here calls dyn_series\sec if in the dyn_series directory
            for k = 1:order
                hpk = dot(vAll(1:k),up(k:-1:1));
                hpfliplr(order-k+1) = hpk;
                vAll(k+1) = 2*dot(hvec(1:k),hpfliplr(order-k+1:order))/k;
                hvec(k+1) = hpk/k;
            end
            h = dyn_series(hvec(1), hvec(2:order+1));
        end
        %=============================
        function h = atan(u)
            %dyn_series/ATAN overloads arctangent of a dyn_series object argument
            order = length(u.coef);
            up = (1:order).*u.coef;
            h0 = atan(u.val);
            v0 = (cos(h0))^(-2);
            hp = zeros(1,order);
            vcoeffliplr = zeros(1,order);
            ufliplr = [fliplr(u.coef) u.val];
            for k = 1:order
                hp(k) = (up(k) - dot(hp(1:k-1),vcoeffliplr(order-(k-2):order)))/v0;
                vpk = 2*dot(up(1:k),ufliplr(order+2-k:order+1));
                vcoeffliplr(order-k+1) = vpk/k;
            end
            h = dyn_series(h0,hp./(1:order));
        end
    end
end