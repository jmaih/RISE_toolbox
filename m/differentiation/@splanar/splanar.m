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
        is_known=true % flag for known functions
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
                    
                    ni=sub_number_of_columns(obj.args{ii});
                    
                    n=max(n,ni);
                    
                end
                
            end
            
            function n=sub_number_of_columns(obj)
                
                n=0;
                
                if isa(obj,'splanar')
                    
                    n=obj.number_of_columns;
                    
                end
                
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
            obj=splanar.prototypize();
            
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
            
            if is_zero(b)|| is_one(a)
                % x^0 = 1
                % 1^x = 1 
	            obj=splanar.prototypize(1);
                
            elseif is_one(b)
                % x^1 = x
                obj=a;
                
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
                
                obj=splanar.prototypize();
                
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
                obj=splanar.prototypize();
                
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
                if strcmp(class(a.args{1}),'splanar')
                    
                    obj=a.args{1};
                    
                else
                    
                    obj=splanar.prototypize(a.args{1});
                    
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
        
        varargout=diff(varargin)
        
        varargout=diff_test(varargin)
        
    end
    
    methods(Hidden)
        
        varargout=intercept_column(varargin)
        
    end
    
    methods(Access=private)
                
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
        
        function [eqtns]=alienize(eqtns,alien_list)
            
            if isempty(alien_list)
                
                return
                
            end
            
            was_char=ischar(eqtns);
            
            if was_char
                
                eqtns=cellstr(eqtns);
                
            end
            
            alien_list=parser.cell2matize(alien_list);
            
            engine=@(x)regexprep(x,['\<',alien_list,'\('],'splanar.backdoor(''$1'',');
            
            for ii=1:numel(eqtns)
                
                eqtn=eqtns{ii};
                
                if isa(eqtn,'function_handle')
                    
                    eqtn=func2str(eqtn);
                    
                    eqtns{ii}=str2func(engine(eqtn));
                    
                else
                    
                    eqtns{ii}=engine(eqtn);
                    
                end
                
            end
            
            if was_char
                
                eqtns=char(eqtns);
                
            end
            
        end
        
        function y=backdoor(f,varargin)
            
            y=splanar(f,varargin);
            
            incid=[];
            
            for ii=1:length(varargin)
                
                arg=varargin{ii};
                
                if isa(arg,'splanar')
                    
                    a_incid=arg.incidence;
                    
                    if isempty(incid)
                    
                    incid=a_incid;
                    
                    elseif ~isempty(a_incid)
                        
                        incid=incid|a_incid;
                        
                    end
                    
                end
                
            end
            
            y.incidence=incid;
            
            y.is_known=false;
            
        end
        
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
        
        varargout=prototypize(varargin)
        
        varargout=rediff(varargin)
        
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

function varargout=splanarize(varargin)

varargout=varargin;

n=nargin;

for iarg=1:n
    
    if ~strcmp(class(varargin{iarg}),'splanar')

	    varargout{iarg}=splanar.prototypize(varargout{iarg});
	    
	end
    
end	

end

function obj=do_trivariate(x,mu,sd,func)

[x,mu,sd]=splanarize(x,mu,sd);

if isnumeric(x.func) && isnumeric(mu.func) && isnumeric(sd.func)
    
    obj=splanar.prototypize(feval(func,x.func,mu.func,sd.func));
    
else
    
    obj=splanar.prototypize(func);
    
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

if isnumeric(a.func) && isnumeric(b.func)
    
    if any(strcmp(func,{'mtimes','mrdivide','mpower'}))
        
        func=func(2:end);
        
    end
    
    obj=splanar.prototypize(feval(func,a.func,b.func));
    
else
    
    obj=splanar.prototypize(func);
    
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

if isnumeric(a.func)
    
    obj=splanar.prototypize(feval(func,a.func));
    
else
    
    obj=splanar.prototypize(func);
    
    obj.args={a};
    
    obj.incidence=a.incidence;
    
end

end
