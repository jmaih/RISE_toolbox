classdef splanar
    % splanar symbolic "planar" differentiation
    %
    % - [abs](splanar/abs)
    % - [acos](splanar/acos)
    % - [acosh](splanar/acosh)
    % - [and](splanar/and)
    % - [asin](splanar/asin)
    % - [asinh](splanar/asinh)
    % - [atan](splanar/atan)
    % - [atanh](splanar/atanh)
    % - [char](splanar/char)
    % - [cos](splanar/cos)
    % - [cosh](splanar/cosh)
    % - [cot](splanar/cot)
    % - [derivatives2functions](splanar/derivatives2functions)
    % - [diff](splanar/diff)
    % - [differentiate](splanar/differentiate)
    % - [eq](splanar/eq)
    % - [erf](splanar/erf)
    % - [exp](splanar/exp)
    % - [ge](splanar/ge)
    % - [get](splanar/get)
    % - [gt](splanar/gt)
    % - [if_elseif](splanar/if_elseif)
    % - [if_then_else](splanar/if_then_else)
    % - [initialize](splanar/initialize)
    % - [intercept_column](splanar/intercept_column)
    % - [is_one](splanar/is_one)
    % - [is_zero](splanar/is_zero)
    % - [isnumeric](splanar/isnumeric)
    % - [kron](splanar/kron)
    % - [le](splanar/le)
    % - [load_varlist](splanar/load_varlist)
    % - [log](splanar/log)
    % - [log10](splanar/log10)
    % - [lt](splanar/lt)
    % - [max](splanar/max)
    % - [min](splanar/min)
    % - [minus](splanar/minus)
    % - [mpower](splanar/mpower)
    % - [mrdivide](splanar/mrdivide)
    % - [mtimes](splanar/mtimes)
    % - [ne](splanar/ne)
    % - [normcdf](splanar/normcdf)
    % - [norminv](splanar/norminv)
    % - [normpdf](splanar/normpdf)
    % - [or](splanar/or)
    % - [plus](splanar/plus)
    % - [power](splanar/power)
    % - [print](splanar/print)
    % - [rdivide](splanar/rdivide)
    % - [set](splanar/set)
    % - [sign](splanar/sign)
    % - [sin](splanar/sin)
    % - [sinh](splanar/sinh)
    % - [splanar](splanar/splanar)
    % - [sqrt](splanar/sqrt)
    % - [tan](splanar/tan)
    % - [tanh](splanar/tanh)
    % - [times](splanar/times)
    % - [uminus](splanar/uminus)
    % - [uplus](splanar/uplus)
    
    properties
        func
        args
        incidence % tells which variable appears
        lineage={} % will serve to trace back all generations of derivatives
        %         prototype
        location % will hold position (column) where the derivative will be stored
    end
    
    properties(Dependent)
        
        number_of_columns=1
        
    end
    
    properties(Constant)
        
        known_functions=do_known_functions()
        
    end
    
    methods
        % constructor
        %------------
        function obj=splanar(f,a)
            % splanar - constructor for splanar objects
            %
            % ::
            %
            %
            %   obj=splanar(f)
            %   obj=splanar(f,a)
            %
            % Args:
            %              %
            %              % - **f** [splanar|char|numeric]: function (string) or variable name or
            %              %   numerical value
            %              %
            %              % - **a** [cell|numeric|splanar]: argument(s) of function f. This is when f
            %              %   is indeed a function
            %              %
            % Returns:
            %    :
            %              %
            %              % - **obj** [splanar]: built object
            %              %
            % Note:
            %              %
            % Example:
            %
            % See also:
            if nargin>0
                
                if strcmp(class(f),'splanar') %#ok<*STISA>
                
                    obj=f;
                
                else
                    
                    obj.func=f;
                    
                    if nargin>1
                    
                        if ~iscell(a)
                        
                            a={a};
                        
                        end
                        
                        obj.args=a;
                    
                    end
                    
                end
                
            end
            
        end
        
        function n=get.number_of_columns(obj)
            
            if isnumeric(obj.func)
                
                n=numel(obj.func);
                
            elseif isempty(obj.args)
                
                n=1;
                
            else
                
                n=1;
                
                for ii=1:numel(obj.args)
                    
                    n=max(n,...
                        sub_number_of_columns(obj.args{ii}));
                    
                end
                
            end
            
            function n=sub_number_of_columns(obj)
                
                n=obj.number_of_columns;
                
            end
            
        end
        
        % overloaded operators
        %---------------------
        function obj=abs(a)
            % abs - overloads abs for splanar
            obj=do_univariate(a,'abs');
            
        end
        
        function obj=acos(a)
            % acos - overloads acos for splanar
            obj=do_univariate(a,'acos');
            
        end
        
        function obj=acosh(a)
            % acosh - overloads acosh for splanar
            obj=do_univariate(a,'acosh');
            
        end
        
        function obj=and(a,b)
            % and - overloads and for splanar
            obj=do_bivariate(a,b,'and');
            
        end
        
        function obj=asin(a)
            % asin - overloads asin for splanar
            obj=do_univariate(a,'asin');
            
        end
        
        function obj=asinh(a)
            % asinh - overloads asinh for splanar
            obj=do_univariate(a,'asinh');
            
        end
        
        function obj=atan(a)
            % atan - overloads atan for splanar
            obj=do_univariate(a,'atan');
            
        end
        
        function obj=atanh(a)
            % atanh - overloads atanh for splanar
            obj=do_univariate(a,'atanh');
            
        end
        
        function obj=betainv(x,a,b)
            % betacdf - overloads betacdf for splanar
            obj=do_trivariate(x,a,b,'betainv');
            
        end
        
        function obj=betacdf(x,a,b)
            % betacdf - overloads betacdf for splanar
            obj=do_trivariate(x,a,b,'betacdf');
            
        end
        
        function obj=betapdf(x,a,b)
            % betapdf - overloads betapdf for splanar
            obj=do_trivariate(x,a,b,'betapdf');
            
        end
        
        function obj=cos(a)
            % cos - overloads cos for splanar
            obj=do_univariate(a,'cos');
            
        end
        
        function obj=cosh(a)
            % cosh - overloads cosh for splanar
            obj=do_univariate(a,'cosh');
            
        end
        
        function obj=cot(a)
            % cot - overloads cot for splanar
            obj=do_univariate(a,'cot');
            
        end
        
        function obj=eq(a,b)
            % eq - overloads eq for splanar
            obj=do_bivariate(a,b,'eq');
            
        end
        
        function obj=erf(a)
            % erf - overloads erf for splanar
            obj=do_univariate(a,'erf');
            
        end
        
        function obj=exp(a)
            % exp - overloads exp for splanar
            obj=do_univariate(a,'exp');
            
        end
        
        function obj=ge(a,b)
            %  - overloads  for splanar
            obj=do_bivariate(a,b,'ge');
            
        end
        
        function obj=gt(a,b)
            % ge - overloads ge for splanar
            obj=do_bivariate(a,b,'gt');
            
        end
        
        function obj=if_elseif(varargin)
            % if_elseif - overloads if_elseif for splanar
            n=nargin;
            
            [varargin{1:n}]=splanarize(varargin{:});
            
            nconst=0;
            
            constants=cell(1,n);
            
            for iarg=1:n
                
                if isnumeric(varargin{iarg})
                    
                    constants{iarg}=varargin{iarg}.func;
                    
                    nconst=nconst+1;
                
                else
                    
                    break
                
                end
                
            end
            % initialize the splanar
            %------------------------
            obj=prototypize(varargin{1});
            
            if nconst==n
                
                obj.func=if_elseif(constants{:});
                
            else
                
                obj.func='if_elseif';
                
                obj.args=varargin;
                
                for iarg=1:n
                    
                    if ~isempty(varargin{iarg}.incidence)
                        
                        if isempty(obj.incidence)
                            
                            obj.incidence=varargin{iarg}.incidence;
                            
                        else
                            
                            obj.incidence=obj.incidence|varargin{iarg}.incidence;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function obj=if_then_else(a,b,c)
            % if_then_else - overloads if_then_else for splanar
            obj=if_elseif(a,b,~a,c);
            
        end
        
        function obj=le(a,b)
            % le - overloads le for splanar
            obj=do_bivariate(a,b,'le');
            
        end
        
        function obj=log(a)
            % log - overloads log for splanar
            obj=do_univariate(a,'log');
            
        end
        
        function obj=log10(a)
            % log10 - overloads log10 for splanar
            obj=log(a)/log(10);
            
        end
        
        function obj=lt(a,b)
            % lt - overloads lt for splanar
            obj=do_bivariate(a,b,'lt');
            
        end
        
        function obj=max(a,b)
            % max - overloads max for splanar
            obj=do_bivariate(a,b,'max');
            
        end
        
        function obj=min(a,b)
            % min - overloads min for splanar
            obj=do_bivariate(a,b,'min');
            
        end
        
        function obj=minus(a,b)
            % minus - overloads minus for splanar
            [a,b]=splanarize(a,b);
            
            if is_zero(b)
                % x - 0 = x
                obj=a;
                
            elseif is_zero(a)
                % 0 - x = - x
                obj=-b;
                
            elseif ischar(b.func) && strcmp(b.func,'uminus')
                % x - (-y) = x + y
                obj=a+b.args{1};
                
            else
                
                obj=do_bivariate(a,b,'minus',false);
                
            end
            
        end
        
        function obj=mpower(a,b)
            % mpower - overloads mpower for splanar
            [a,b]=splanarize(a,b);
            
            obj=prototypize(a);
            
            if is_zero(b)
                % x^0 = 1
                obj.func=1;
                
            elseif is_one(b)
                % x^1 = x
                obj=a;
                
            elseif is_one(a)
                % 1^x = 1
                obj.func=1;
                
            else
                
                obj=do_bivariate(a,b,'mpower',false);
                
            end
            
        end
        
        function obj=mrdivide(a,b)
            % mrdivide - overloads mrdivide for splanar
            [a,b]=splanarize(a,b);
            
            if any(b.func==0)
                
                % x/0 impossible
                error('dividing by zero not allowed')
                
            elseif is_one(b)
                % x/1 = x
                obj=a;
                
            elseif is_zero(a)
                
                obj=prototypize(a);
                
                obj.func=0;
                
            elseif strcmp(a.func,'uminus')
                
                if strcmp(b.func,'uminus')
                    % (-x) / (-y) = x * y
                    obj=mrdivide(a.args{1},b.args{1});
                    
                else
                    % (-x) / y = -(x*y)
                    obj=-mrdivide(a.args{1},b);
                    
                end
                
            else
                
                obj=do_bivariate(a,b,'mrdivide',false);
                
            end
            
        end
        
        function obj=mtimes(a,b)
            
            % mtimes - overloads mtimes for splanar
            [a,b]=splanarize(a,b);
            
            isnum_a=isnumeric(a);
            
            isnum_b=isnumeric(b);
            
            if strcmp(a.func,'uminus')
                
                if strcmp(b.func,'uminus')
                    % (-x) * (-y) = x * y
                    obj=mtimes(a.args{1},b.args{1});
                    
                else
                    
                    % (-x) * y = -(x*y)
                    obj=-mtimes(a.args{1},b);
                    
                end
                
            elseif strcmp(b.func,'uminus')
                
                % x * (-y) = -(x*y)
                obj=-mtimes(a,b.args{1});
                
            elseif (isnum_a && is_zero(a))||(isnum_b && is_zero(b))
                % 0 * x = 0
                % x * 0 = 0
                obj=prototypize(a);
                
                obj.func=0;% obj.func=0*ones(1,maxcols);
                
            elseif isnum_a && is_one(a) % maxcols==1 &&
                % 1 * x = x
                obj=b;
                
            elseif isnum_b && is_one(b)% maxcols==1 &&
                % x * 1 = x
                obj=a;
                
            else
                
                obj=do_bivariate(a,b,'mtimes',false);
                
            end
            
        end
        
        function obj=ne(a,b)
            % ne - overloads ne for splanar
            obj=do_bivariate(a,b,'ne');
            
        end
        
        function obj=norminv(x,mu,sd)
            % normcdf - overloads normcdf for splanar
            if nargin<3
                
                sd=1;
                
                if nargin<2
                    
                    mu=0;
                    
                end
                
            end
            
            obj=do_trivariate(x,mu,sd,'norminv');
            
        end
        
        function obj=normcdf(x,mu,sd)
            % normcdf - overloads normcdf for splanar
            if nargin<3
                
                sd=1;
                
                if nargin<2
                    
                    mu=0;
                    
                end
                
            end
            
            obj=do_trivariate(x,mu,sd,'normcdf');
            
        end
        
        function obj=normpdf(x,mu,sd)
            % normpdf - overloads normpdf for splanar
            if nargin<3
                
                sd=1;
                
                if nargin<2
                    
                    mu=0;
                    
                end
                
            end
            
            obj=do_trivariate(x,mu,sd,'normpdf');
            
        end
        
        function obj=or(a,b)
            % or - overloads or for splanar
            obj=do_bivariate(a,b,'or');
            
        end
        
        function obj=plus(a,b)
            % plus - overloads plus for splanar
            [a,b]=splanarize(a,b);
            
            if ischar(b.func) && strcmp(b.func,'uminus')
                % x + (-y) = x - y
                obj=minus(a,b.args{1});
                
            elseif is_zero(a)
                % 0 + x = x
                obj=b;
                
            elseif is_zero(b)
                % x + 0 = x
                obj=a;
                
            elseif isequal(a,b) %% this is going to be expensive ... a.pointer==b.pointer
                % x + x = 2x
                obj=2*a;
                
            else
                
                obj=do_bivariate(a,b,'plus',false);
                
            end
            
        end
        
        function obj=sign(a)
            % sign - overloads sign for splanar
            obj=do_univariate(a,'sign');
            
        end
        
        function obj=sin(a)
            % sin - overloads sin for splanar
            obj=do_univariate(a,'sin');
            
        end
        
        function obj=sinh(a)
            % sinh - overloads sinh for splanar
            obj=do_univariate(a,'sinh');
            
        end
        
        function obj=sqrt(a)
            % sqrt - overloads sqrt for splanar
            obj=do_univariate(a,'sqrt');
            
        end
        
        function obj=steady_state(a)
            % sqrt - overloads sqrt for splanar
            obj=do_univariate(a,'steady_state');
            
        end
        
        function obj=tan(a)
            % tan - overloads tan for splanar
            obj=do_univariate(a,'tan');
            
        end
        
        function obj=tanh(a)
            % tanh - overloads tanh for splanar
            obj=do_univariate(a,'tanh');
            
        end
        
        function obj=uminus(a)
            % uminus - overloads uminus for splanar
            if strcmp(a.func,'uminus')
                
                % Simplify -(-x) in x
                obj=prototypize(a);
                
                if strcmp(class(a.args{1}),'splanar')
                    
                    obj=a.args{1};
                    
                else
                    
                    obj.func=a.args{1};
                    
                end
                
            else
                
                obj=do_univariate(a,'uminus');
                
            end
            
        end
        
        function obj=uplus(a)
            % uplus - overloads uplus for splanar
            obj=a;
            
        end
        % bizarre functions
        %------------------
        function obj=power(varargin)
            % power - overloads power for splanar
            obj=mpower(varargin{:});
            
        end
        
        function obj=rdivide(varargin)
            % rdivide - overloads rdivide for splanar
            obj=mrdivide(varargin{:});
            
        end
        
        function obj=times(varargin)
            % times - overloads times for splanar
            obj=mtimes(varargin{:});
            
        end
        
        % utility functions
        %------------------
        varargout=char(varargin)
        
        varargout=kron(varargin)
        
        varargout=set(varargin)
        
        varargout=get(varargin)
        
        function d=diff(x,wrt,pointer)
            % diff - overloads diff for splanar
            if nargin<3
                
                pointer=[];
                
            end
            
            
            [nrows,ncols]=size(x);
            
            if nrows>1||ncols>1
                
                d=cell(nrows,ncols);
                
                for irow=1:nrows
                    
                    for icol=1:ncols
                        
                        d{irow,icol}=diff(x(irow,icol),wrt);
                        
                    end
                    
                end
                
                return
                
            end
            
            if isempty(x.incidence)
                % numbers/vectors and variables which are not part of differentiation
                % automatically receive 0 as derivative
                d=prototypize(x);d.func=0;
                
            elseif isempty(x.args)
                % variables that are part of differentiation
                d=prototypize(x);
                
                d.func=double(x.incidence(wrt));
                % why do we need to double it?
                % we need to carry around vectors in order to be able to do
                % higher-order derivatives
            else
                % functions of variables
                if_elseif_flag=strcmp(x.func,'if_elseif');
                
                if_then_else_flag=strcmp(x.func,'if_then_else');
                
                nargs=numel(x.args);
                
                d_args=x.args;
                
                for iarg=1:nargs
                    
                    if ~isempty(pointer) && isnumeric(x.args{iarg}) && numel(x.args{iarg}.func)>1
                        
                        x.args{iarg}.func=x.args{iarg}.func(pointer);
                    
                    end
                    
                    if (if_elseif_flag && rem(iarg,2))||(iarg==1 && if_then_else_flag)
                        
                        continue
                        
                    end
                    
                    d_args{iarg}=diff(x.args{iarg},wrt,pointer);
                    
                end
                % compose derivatives
                %--------------------
                switch x.func
                    
                    case 'plus'
                        
                        d=d_args{1}+d_args{2};
                        
                    case 'minus'
                        
                        d=d_args{1}-d_args{2};
                        
                    case 'mtimes'
                        
                        t11 = d_args{1}*x.args{2};
                        
                        t12 = d_args{2}*x.args{1};
                        
                        d=t11+t12;
                        
                    case 'mrdivide'
                        
                        if ~is_zero(d_args{2})
                            
                            t11 = d_args{1}*x.args{2};
                            
                            t12 = d_args{2}*x.args{1};
                            
                            t13 = t11-t12;
                            
                            t14 = x.args{2}^2;
                            
                            d=t13/t14;
                            
                        else
                            
                            d=d_args{1}/x.args{2};
                            
                        end
                        
                    case {'lt','gt','le','ge','eq','ne','or','and','sign',...
                            'steady_state'}
                        
                        % derivatives are zero
                        d=prototypize(x); d.func=0;
                        
                    case 'if_then_else'
                        
                        d=if_then_else(x.args{1},d_args{2},d_args{3});
                        
                    case 'if_elseif'
                        
                        the_args=x.args;
                        
                        the_args(2:2:end)=d_args(2:2:end);
                        
                        d=if_elseif(the_args{:});
                        
                    case 'mpower'
                        
                        t11 = log(x.args{1});
                        
                        t12 = d_args{2}*t11;
                        
                        t13 = d_args{1}*x.args{2};
                        
                        d= t12*x+t13*x.args{1}^(x.args{2}-1);
                        %         t14 = t13/x.args{1}; % <-- creates problems when x.args{1}=0
                        %         t15 = t12+t14;
                        %         d= t15*x.args{1}^x.args{2};
                    case 'max'
                        
                        t11 = x.args{1}>x.args{2};
                        
                        t12 = t11*d_args{1};
                        
                        t13 = 1-t11;
                        
                        t14 = t13*d_args{2};
                        
                        d= t14+t12;
                        
                    case 'min'
                        
                        t11 = x.args{2}>x.args{1};
                        
                        t12 = t11*d_args{1};
                        
                        t13 = 1-t11;
                        
                        t14 = t13*d_args{2};
                        
                        d= t14+t12;
                        
                    case 'uminus'
                        
                        d= -d_args{1};
                        
                    case 'uplus'
                        
                        d= d_args{1};
                        
                    case 'exp'
                        
                        d= d_args{1}*x;
                        
                    case 'log'
                        
                        d= d_args{1}/x.args{1};
                        
                    case 'log10'
                        
                        t11 = exp(1);
                        
                        t12 = log10(t11);
                        
                        t13 = d_args{1}/x.args{1};
                        
                        d= t12*t13;
                        
                    case 'cos'
                        
                        t11 = sin(x.args{1});
                        
                        t12 = -t11;
                        
                        d= d_args{1}*t12;
                        
                    case 'sin'
                        
                        t11 = cos(x.args{1});
                        
                        d= d_args{1}*t11;
                        
                    case 'tan'
                        
                        d= d_args{1}*(1+x^2);
                        
                    case 'cot'
                        
                        d= -d_args{1}*(1-x^2); % <--- d=-d_args{1}/sin(x.args{1})^2;
                        
                    case 'acos'
                        
                        d= -d_args{1}/sqrt(1-x^2);
                        
                    case 'asin'
                        
                        t11 = cos(x);
                        
                        d= d_args{1}/t11;
                        
                    case 'atan'
                        
                        t11 = x.args{1}^2;
                        
                        t12 = 1+t11;
                        
                        d= d_args{1}/t12;
                        
                    case 'cosh'
                        
                        t11 = sinh(x.args{1});
                        
                        d= d_args{1}*t11;
                        
                    case 'sinh'
                        
                        t11 = cosh(x.args{1});
                        
                        d= d_args{1}*t11;
                        
                    case 'tanh'
                        
                        d= d_args{1}*(1-x^2);
                        
                    case 'acosh'
                        
                        t11 = sinh(x);
                        
                        d= d_args{1}/t11;
                        
                    case 'asinh'
                        
                        t11 = cosh(x);
                        
                        d= d_args{1}/t11;
                        
                    case 'atanh'
                        
                        t11 = x.args{1}^2;
                        
                        t12 = 1-t11;
                        
                        d= d_args{1}*t12;
                        
                    case 'sqrt'
                        
                        t11 = 2*x;
                        
                        d= d_args{1}/t11;
                        
                    case 'abs'
                        
                        t11 = sign(x.args{1});
                        
                        d= t11*d_args{1};
                        
                    case 'erf'
                        % x^2
                        t11 = mpower(x.args{1},2);
                        
                        % exp(x^2)
                        t12 =  exp(t11);
                        
                        % sqrt(pi)
                        t11 = sqrt(pi);
                        
                        % sqrt(pi)*exp(x^2)
                        t13 = t11*t12;
                        
                        % 2/(sqrt(pi)*exp(x^2));
                        t14 = 2/t13;
                        
                        % (2/(sqrt(pi)*exp(x^2)))*dx;
                        d= t14*d_args{1};
                        
                    case 'norminv'
                        
                        d=(d_args{1}+d_args{2}+d_args{3})/normpdf(...
                            norminv(x.args{1},x.args{2},x.args{3}));
                        
                    case 'normcdf'
                        
                        d=(d_args{1}+d_args{2}+d_args{3})*normpdf(x.args{1},x.args{2},x.args{3});
                    
                    case 'normpdf'
                        
                        y=(x.args{1}-x.args{2})/x.args{3};
                        
                        ss= d_args{3}/x.args{3};
                        
                        d=((ss*y-(d_args{1}-d_args{2})/x.args{3})*y-ss)*x;
                        
                    case 'betainv'
                        
                        d=(d_args{1}+d_args{2}+d_args{3})/betapdf(...
                            betainv(x.args{1},x.args{2},x.args{3}));
                        
                    case 'betacdf'
                        
                        d=(d_args{1}+d_args{2}+d_args{3})*betapdf(x.args{1},x.args{2},x.args{3});
                   
                    case 'betapdf'
                        % https://en.wikipedia.org/wiki/Beta_distribution
                        % this assumes that the derivatives wrt a and b are
                        % zero !
                        numerator=(x.args{2}+x.args{3}-2)*x.args{1};
                        
                        denominator= (x.args{1}-1)*x.args{1};
                        
                        d=numerator/denominator*betapdf(x.args{1},x.args{2},x.args{3});
                
                end
                
            end
            
        end
        
    end
    
    methods(Access=private)
        
        function obj=intercept_column(obj,pointer)
            % intercept_column - builds a scalar splanar object from a
            % vectorized splanar object. The pointer argument points the
            % element in the vector to be used.
            if obj.number_of_columns>1
                
                if isnumeric(obj) && numel(obj.func)>1
                    
                    obj.func=obj.func(pointer);
                    
                elseif ~isempty(obj.args)
                    
                    for iarg=1:numel(obj.args)
                        
                        if ~isa(obj.args{iarg},'splanar')
                            
                            continue
                            
                        end
                        
                        obj.args{iarg}=intercept_column(obj.args{iarg},pointer);
                        
                    end
                    
                    % the line below will correct for zeros, ones, etc. as well as
                    % incidences.
                    obj=feval(obj.func,obj.args{:});
                    
                end
                
            end
            
        end
        
        varargout=load_varlist(varargin)
        
        function flag=isnumeric(a)
            
            % isnumeric - checks whether a splanar object is numeric
            flag=isnumeric(a.func)||islogical(a.func);
            
        end
        
        function flag=is_zero(a)
            
            % is_zero - checks whether a splanar object is zero
            afunc=a.func;
            
            flag=(isnumeric(afunc)||islogical(afunc)) && all(afunc==0);
            
        end
        
        function flag=is_one(a)
            % is_one - checks whether a splanar object is 1
            
            afunc=a.func;
            
            flag=(isnumeric(afunc)||islogical(afunc)) && all(afunc==1);
            
        end
        
    end
    
    methods(Static)
        
        function var_list=initialize(var_list,wrt_list)
            % initialize - initializes splanar objects for differentiation
            if ischar(var_list)
                
                var_list=cellstr(var_list);
                
            end
            
            if ischar(wrt_list)
                
                wrt_list=cellstr(wrt_list);
                
            end
            
            nwrt=numel(wrt_list);
            
            nvlist=numel(var_list);
            
            if numel(unique(wrt_list))~=nwrt
                
                error('repeated variable names in wrt_list')
                
            end
            
            if numel(unique(var_list))~=nvlist
                
                error('repeated variable names in var_list')
                
            end
            
            proto_incidence=false(1,nwrt);
            %             proto_=splanar();
            for ivar=1:nvlist
                
                vname=var_list{ivar};
                
                loc=find(strcmp(vname,wrt_list));
                
                var_list{ivar}=splanar(vname);
                
                if ~isempty(loc)
                    
                    incid_=proto_incidence;
                    
                    incid_(loc)=true;
                    
                    var_list{ivar}=set(var_list{ivar},'incidence',sparse(incid_));
                
                end
                %                 var_list{ivar}=set(var_list{ivar},'prototype',proto_);
            end
            
        end
        
        varargout=differentiate(varargin)
        
        varargout=print(varargin)
        
    end
    
end

function kf=do_known_functions()

kf={'abs'
    'acos'
    'acosh'
    'and'
    'asin'
    'asinh'
    'atan'
    'atanh'
    'betacdf'
    'betainv'
    'betapdf'
    'char'
    'cos'
    'cosh'
    'cot'
%     'diff'
%     'differentiate'
    'eq'
    'erf'
    'exp'
    'ge'
    'get'
    'gt'
    'if_elseif'
    'if_then_else'
    'initialize'
    'kron'
    'le'
    'log'
    'log10'
    'lt'
    'max'
    'min'
    'minus'
    'mpower'
    'mrdivide'
    'mtimes'
    'ne'
    'normcdf'
    'norminv'
    'normpdf'
    'or'
    'plus'
    'power'
    'print'
    'rdivide'
    'set'
    'sign'
    'sin'
    'sinh'
    'splanar'
    'sqrt'
    'steady_state'
    'tan'
    'tanh'
    'times'
    'uminus'
    'uplus'
}.';

end

function obj=prototypize(obj)

obj.func=[];

obj.args=[];

obj.incidence=[];

end

function varargout=splanarize(varargin)

varargout=varargin;

n=nargin;

obj=[];

guy_is_planar=false(1,n);

for iarg=1:n
    
    guy_is_planar(iarg)=strcmp(class(varargin{iarg}),'splanar');
    
    if isempty(obj) && guy_is_planar(iarg)
        
        obj=prototypize(varargin{iarg});
        
    end
    
end

for iarg=find(~guy_is_planar)
    
    tmp=varargout{iarg};
    
    varargout{iarg}=obj;
    
    varargout{iarg}.func=tmp;
    
end

end

function obj=do_trivariate(x,mu,sd,func)
[x,mu,sd]=splanarize(x,mu,sd);
% initialize the splanar
%------------------------
obj=prototypize(x);

if isnumeric(x.func) && isnumeric(mu.func) && isnumeric(sd.func)
    
    obj.func=feval(func,x.func,mu.func,sd.func);
    
else
    
    obj.func=func;
    
    obj.args={x,mu,sd};
    
    obj.incidence=x.incidence;
    
    if ~isempty(mu.incidence)
        
        if isempty(obj.incidence)
            
            obj.incidence=mu.incidence;
            
        else
            
            obj.incidence=obj.incidence|mu.incidence;
            
        end
        
    end
    
    if ~isempty(sd.incidence)
        
        if isempty(obj.incidence)
            
            obj.incidence=sd.incidence;
            
        else
            
            obj.incidence=obj.incidence|sd.incidence;
            
        end
        
    end
    
end

end

function obj=do_bivariate(a,b,func,re_splanarize)

if nargin<4
    
    re_splanarize=true;
    
end

if re_splanarize
    
    [a,b]=splanarize(a,b);
    
end
% initialize the splanar
%------------------------
obj=prototypize(a);

if isnumeric(a.func) && isnumeric(b.func)
    
    if any(strcmp(func,{'mtimes','mrdivide','mpower'}))
        
        func=func(2:end);
        
    end
    
    obj.func=feval(func,a.func,b.func);
    
else
    
    obj.func=func;
    
    obj.args={a,b};
    
    obj.incidence=a.incidence;
    
    if ~isempty(b.incidence)
        
        if isempty(obj.incidence)
            
            obj.incidence=b.incidence;
            
        else
            
            obj.incidence=obj.incidence|b.incidence;
            
        end
        
    end
    
end

end

function obj=do_univariate(a,func)
% initialize the splanar
%------------------------
obj=prototypize(a);

if isnumeric(a.func)
    
    obj.func=feval(func,a.func);
    
else
    
    obj.func=func;
    
    obj.args={a};
    
    obj.incidence=a.incidence;
    
end

end
