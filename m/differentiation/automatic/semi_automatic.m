classdef semi_automatic
    properties
        x
        dx
    end
    methods
        function obj=semi_automatic(a,b)
            if nargin
                obj.x=a;
                if nargin>1
                    obj.dx=b;
                    if nargin>2
                        error([mfilename,':: number of arguments cannot exceed 2'])
                    end
                else
                    obj.dx=0;
                end
            end
        end
%         function [x,dx] = semi_automatic_get(ad,mode)
%             % ADIFFGET returns parts of an semi_automatic object.
%             %
%             % x      = semi_automaticget(A,'value')      returns the value of the semi_automatic object A
%             % dx     = semi_automaticget(A,'derivative') returns the jacobian
%             % [x,dx] = semi_automaticget(A)              returns both the value and the jacobian
%             %
%             % The value of an semi_automatic object is what value the object would have had it just
%             % been a vector. The jacobian is the derivative of the value with respect to the
%             % variables implicitly created when the semi_automatic object was first created.
%             %
%             % For example, a=semi_automatic(1:10) creates an semi_automatic object of 10 variables, with values 1 to 10.
%             % [x,dx] = semi_automaticget(a) will then yield x=(1:10)', and dx is a 10 by 10 sparse identity matrix.
%             % If we then compute b = sum(a), b is an semi_automatic object (because sum is defined for them).
%             % [x,dx] = semi_automaticget(b) will yield x=55 i.e. sum(1:10), and dx =ones(1,10), since that is the
%             % jacobian matrix for the sum.
%             
%             if nargin==1
%                 x = ad.x;
%                 dx = ad.dx;
%             else
%                 switch mode(1)
%                     case 'v'
%                         x = ad.x;
%                     case 'd'
%                         x = ad.dx;
%                 end
%             end
%         end
        function c = mldivide(a,b)
            % MLDIVIDE implements a\b, where a is double and b is an semi_automatic object, or
            % when a is an semi_automatic object with length 1.
            
            switch [class(a),class(b)]
                case 'doublesemi_automatic'
                    c = semi_automatic(a\b.x, a\b.dx);
                case {'semi_automaticdouble','semi_automaticsemi_automatic'}
                    if length(a)==1
                        c = times(b,1./a);
                    else
                        error('Cannot do a\b when a is an semi_automatic object with length>1')
                    end
                otherwise
                    error(['\ not defined for ',class(a),' and ',class(b)]);
            end
        end
        function c = minus(a,b)
            % MINUS implements a-b, where either a or b is an semi_automatic object.
            
            c = plus(a,-b);
        end
        function c = min(a,b)
            % MIN for semi_automatic objects. This selects the value of a which is a maximum.
            % Both min(a) and min(a,b) will work. min(a,[],dim) is meaningless because
            % the semi_automatic object is a column vector; use min(a) instead.
            
            if nargin==1
                [y,i] = min(a.x);
                c = semi_automatic(y,a.dx(i,:));
            else
                switch [class(a),class(b)]
                    case 'doublesemi_automatic'
                        ok = b.x<a;
                        c = semi_automatic(min(a,b.x), (ok*b.dx));
                    case 'semi_automaticdouble'
                        ok = a.x<b;
                        c = semi_automatic(min(a.x,b), (ok*a.dx));
                    case 'semi_automaticsemi_automatic'
                        if size(a.dx,1)~=size(b.dx,1)
                            if size(a.dx,1)==1
                                a.dx = repmat(a.dx,size(b.dx,1),1);
                            elseif size(b.dx,1)==1
                                b.dx = repmat(b.dx,size(a.dx,1),1);
                            end
                        end
                        ok = a.x<b.x;
                        c = semi_automatic(min(a.x,b.x), (ok*a.dx)+(1-ok*b.dx));
                    otherwise
                        error(['Can''t compute max of ',class(a),' and ',class(b)]);
                end
            end
        end
        function c = max(a,b)
            % MAX for semi_automatic objects. This selects the value of a which is a maximum.
            % Both max(a) and max(a,b) will work. max(a,[],dim) is meaningless because
            % the semi_automatic object is a column vector; use max(a) instead.
            
            if nargin==1
                [y,i] = max(a.x);
                c = semi_automatic(y,a.dx(i,:));
            else
                switch [class(a),class(b)]
                    case 'doublesemi_automatic'
                        ok = b.x>a;
                        c = semi_automatic(max(a,b.x), (ok*b.dx));
                    case 'semi_automaticdouble'
                        ok = a.x>b;
                        c = semi_automatic(max(a.x,b), (ok*a.dx));
                    case 'semi_automaticsemi_automatic'
                        if size(a.dx,1)~=size(b.dx,1)
                            if size(a.dx,1)==1
                                a.dx = repmat(a.dx,size(b.dx,1),1);
                            elseif size(b.dx,1)==1
                                b.dx = repmat(b.dx,size(a.dx,1),1);
                            end
                        end
                        ok = a.x>b.x;
                        c = semi_automatic(max(a.x,b.x), (ok*a.dx)+(1-ok*b.dx));
                    otherwise
                        error(['Can''t compute max of ',class(a),' and ',class(b)]);
                end
            end
        end
        function c = log10(a)
            % LOG10 for semi_automatic objects.
            
            c = semi_automatic(log(a.x)/log(10), (1./a.x/log(10)*a.dx));
        end
        function c = log(a)
            % LOG for semi_automatic objects.
            
            c = semi_automatic(log(a.x),(1./a.x*a.dx));
        end
        function c = isreal(a)
            % ISREAL for semi_automatic objects
            c = isreal(a.x);
        end
        function c = exp(a)
            % EXP for semi_automatic objects.
            
            c = semi_automatic(exp(a.x),(exp(a.x)*a.dx));
        end
        function d = double(a)
            d = a.x;
        end
        function da = diff(ad,n)
            % DIFF for semi_automatic objects, returns the discrete derivative
            % along the value of the semi_automatic object.
            
            if nargin<2, n=1; end
            da = semi_automatic(diff(ad.x,n),diff(ad.dx,n));
        end
        function c = ctranspose(ad)
            % CTRANSPOSE returns the transpose of an semi_automatic object.
            % Unlike vectors, the transpose must be used with care, and
            % is really only safe in places like bilinear or quadratic forms,
            % e.g. a'*B*c
            
            c = semi_automatic(ad.x',ad.dx);
        end
        function c = cosh(a)
            % COSH of an semi_automatic object.
            
            c = semi_automatic(cosh(a.x),(sinh(a.x)*a.dx));
            
        end
        function c = cos(a)
            % COS of an semi_automatic object.
            
            c = semi_automatic(cos(a.x),(-sin(a.x)*a.dx));
            
        end
        function c = conj(a)
            % CONJ for semi_automatic objects
            c=semi_automatic(conj(a.x),conj(a.dx));
        end
        function c = atan(a)
            % ATAN for semi_automatic objects.
            
            c = semi_automatic(atan(a.x),(1./sqrt(1+a.x.^2)*a.dx));
        end
        function c = asin(a)
            % ASIN for semi_automatic objects.
            
            c = semi_automatic(asin(a.x),(1./sqrt(1-a.x.^2)*a.dx));
        end
        function c = acos(a)
            % ACOS for semi_automatic objects.
            c = semi_automatic(acos(a.x),(-1./sqrt(1-a.x.^2)*a.dx));
        end
        function c = abs(ad)
            % ABS for semi_automatic objects
            
            if all(isreal(ad.x))&& all(isreal(ad.dx(:)))
                c = semi_automatic(abs(ad.x),(sign(ad.x)*ad.dx));
            else
                c = sqrt(ad.*conj(ad));
            end
            
        end
        function c = uminus(ad)
            % UMINUS implements unary minus for semi_automatic objects
            
            % CHECKED
            
            c = semi_automatic( -ad.x, -ad.dx);
        end
        function c = times(a,b)
            % TIMES implements a.*b, where either a or b is an semi_automatic object.
            
            switch [class(a),class(b)]
                
                case 'semi_automaticdouble'
                    c = semi_automatic(b.*a.x, (b*a.dx));
                    
                case 'doublesemi_automatic'
                    c = semi_automatic(a.*b.x, (a*b.dx));
                    
                case 'semi_automaticsemi_automatic'
                    c = semi_automatic( a.x.*b.x, b.x*a.dx+a.x*b.dx);
                    
                otherwise
                    error(['Can''t multiply (.*) ',class(a),' and ',class(b)]);
                    
            end
            
        end
        function c = tanh(a)
            % TANH for semi_automatic objects.
            
            c = semi_automatic( tanh(a.x), 1./cosh(a.x).^2*a.dx);
        end
        function c = tan(a)
            % TAN for semi_automatic objects.
            
            c = semi_automatic(tan(a.x), 1./cos(a.x).^2*a.dx);
            
        end
        function c = sum(a,b)
            % SUM of the semi_automatic object.
            if nargin==1
            c = semi_automatic(sum([a.x]), sum([a.dx]));
            else
                c=plus(a,b);
            end
        end
        function c = sqrt(a)
            % SQRT for semi_automatic objects.
            
            c = a.^(0.5);
        end
        function c = sinh(a)
            % SINH for semi_automatic objects.
            
            c = semi_automatic( sinh(a.x), (cosh(a.x)*a.dx));
        end
        function c = real(a)
            % REAL for semi_automatic objects
            c=semi_automatic(real(a.x),real(a.dx));
        end
        function c = rdivide(a,b)
            % RDIVIDE implements a./b, where either a or b is an semi_automatic object.
            switch [class(a),class(b)]
                case 'semi_automaticdouble'
                    c = times(a,1./b);
                    
                case {'doublesemi_automatic','semi_automaticsemi_automatic'}
                    c = times(a,reciprocal(b));
                    
                otherwise
                    error(['Can''t divide (./) ',class(a),' and ',class(b)]);
            end
            
        end
        function c = power(a,b)
            % POWER implements a.^b, where either a or bor both is an semi_automatic object.
            
            switch [class(a),class(b)]
                
                case 'semi_automaticdouble'
                    c = semi_automatic( a.x.^b, (b.*a.x.^(b-1)*a.dx));
                case 'doublesemi_automatic'
                    c = semi_automatic( a.^b.x, a.^b.x.*log(a)*b.dx);
                case 'semi_automaticsemi_automatic'
                    c = exp(log(a).*b);
                    % Commented out code does it the hard way:
                    % cx  = a.x.^b.x;
                    % cdx = (cx.*log(a.x), b.dx) + (a.x.^(b.x-1), a.dx);
                    % c   = class(struct('x',cx,'dx',cdx),'semi_automatic');
                    
                otherwise
                    error(['Can''t do ',class(a),'.^',class(b)]);
                    
            end
            
        end
        function c = plus(a,b)
            % PLUS implements a+b, where either a or b is an semi_automatic object.
            
            switch [class(a),class(b)]
                
                case 'semi_automaticdouble'
                    c = semi_automatic(a.x+b, a.dx);
                    
                case 'doublesemi_automatic'
                    c = semi_automatic(b.x+a, b.dx);
                    
                case 'semi_automaticsemi_automatic'
                    c = semi_automatic(a.x+b.x,a.dx+b.dx);
                    
                otherwise
                    error(['Can''t add ',class(a),' and ',class(b)]);
            end
        end
        function n = norm(a)
            n = sqrt(sum(a.^2));
        end
        function c = mtimes(a,b)
            % MTIMES implements a*b, where one of a or b is an semi_automatic object.
            % If both are semi_automatic objects, this will only work if at least one
            % has unit length, or a is a transposed semi_automatic object. Otherwise there
            % will be an error of some sort.
            
            if all(size(a)==1) || all(size(b)==1)
                c = times(a,b);
                return
            end
            
            switch [class(a),class(b)]
                
                case 'semi_automaticdouble'
                    c = semi_automatic(a.x*b,b'*a.dx);
                    
                case 'doublesemi_automatic'
                    c = semi_automatic(a*b.x,a*b.dx);
                    
                case 'semi_automaticsemi_automatic'
                    % this will only work if a was transposed
                    c = semi_automatic(a.x*b.x, a.x*b.dx+b.x'*a.dx);
                    
                otherwise
                    error(['Can''t matrix multiply ',class(a),' and ',class(b)]);
                    
            end
        end
        function c = mrdivide(a,b)
            % MRDIVIDE implements a/b, where either a or b is an semi_automatic object.
            % It is mapped to a./b where possible (i.e. when b has length 1), and
            % yields an error otherwise.
            
            if numel(b)==1
                c = times(a,1./b);
            else
                error('Cannot use / for semi_automatic objects when denominator is nonscalar');
            end
        end
        function c = mpower(a,b)
            % MPOWER implements c = a^b, where either a or b is an semi_automatic object.
            % For semi_automatic objects, a^b is taken to be equivalent to a.^b
            
            c=power(a,b);
        end
        
    end
    methods(Access=private)
        function ad1 = reciprocal(ad)
            % RECIPROCAL computes the reciprocal of an adiff object,
            % used for division.
            
            ad1 = semi_automatic(1./ad.x,-1./ad.x.^2*ad.dx);
        end
    end
end

