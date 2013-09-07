classdef planar
%     varlist={'a','b','c','x','d','e'}; wrt={'b','x','e'};
%     wrt_id=1:numel(wrt);
%     vlist=planar.initialize(varlist,wrt);
%     testfunc=@(a,b,c,x,d,e)exp(a*x^2+b*x+c+log(d/sin(e)));
%     pp=testfunc(vlist{:});
%     d=diff(pp,wrt_id);
%     c=char(d)
%     cc=char(d,1)

    properties
        func
        args
    end
    properties(SetAccess=protected)
        incidence
    end
    methods
        function obj=planar(f,a)
            % symbolic differentiation in vectorized form.
            % f can be:
            %         - a planar object
            %         - the name of an overloaded function
            %         - the name of a variable
            %         - a scalar or a vector of numbers
            % a, if supplied, is a cell containing the arguments of the
            % function. Those arguments must be planar objects
            if isa(f,'planar')
                obj=f;
                return
            end
            if nargin>0
                obj.func=f;
                if nargin>1
                    if ~iscell(a)
                        a={a};
                    end
                    obj.args=a;
                end
            end
        end
        varargout=abs(varargin)
        varargout=acos(varargin)
        varargout=acosh(varargin)
        varargout=and(varargin)
        varargout=asin(varargin)
        varargout=asinh(varargin)
        varargout=atan(varargin)
        varargout=atanh(varargin)
        varargout=char(varargin)
        varargout=cos(varargin)
        varargout=cosh(varargin)
        varargout=cot(varargin)
        varargout=diff(varargin)
        varargout=eq(varargin)
        varargout=erf(varargin)
        varargout=exp(varargin)
        varargout=ge(varargin)
        varargout=gt(varargin)
        varargout=if_elseif(varargin)
        varargout=if_then_else(varargin)
        varargout=is_one(varargin)
        varargout=is_zero(varargin)
        varargout=isnumeric(varargin)
        varargout=kron(varargin)
        varargout=le(varargin)
        varargout=load_varlist(varargin)
        varargout=log(varargin)
        varargout=log10(varargin)
        varargout=lt(varargin)
        varargout=max(varargin)
        varargout=min(varargin)
        varargout=minus(varargin)
        varargout=mpower(varargin)
        function varargout=power(varargin)
            varargout=cell(1,nargout);
            [varargout{1:nargout}]=mpower(varargin{:});
        end
        varargout=mrdivide(varargin)
        function varargout=rdivide(varargin)
            varargout=cell(1,nargout);
            [varargout{1:nargout}]=mrdivide(varargin{:});
        end
        varargout=mtimes(varargin)
        function varargout=times(varargin)
            varargout=cell(1,nargout);
            [varargout{1:nargout}]=mtimes(varargin{:});
        end
        varargout=ne(varargin)
        varargout=normalcdf(varargin)
        varargout=normalpdf(varargin)
        varargout=or(varargin)
        varargout=plus(varargin)
        varargout=recompose(varargin)
        varargout=set(varargin)
        varargout=sign(varargin)
        varargout=sin(varargin)
        varargout=sinh(varargin)
        varargout=tan(varargin)
        varargout=uminus(varargin)
        varargout=uplus(varargin)
    end
    methods(Static)
        varargout=initialize(varargin)
        varargout=multinary_operation(varargin)
        varargout=commute(varargin)
    end
end