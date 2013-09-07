classdef rise_sym < handle
    properties
        func
        args
        incidence
        ref
        key
        ncalls=0
    end
    methods
        function obj=rise_sym(func,args,incidence)
            if nargin>0
                if strcmp(class(func),'rise_sym')%#ok<STISA> %isa(func,'rise_sym')
                    obj=func;
                else
                    obj.func=func;
                    if nargin>1
                        if ~iscell(args)
                            args={args};
                        end
                        obj.args=args;
                        if nargin>2
                            obj.incidence=incidence;
                        end
                    end
                end
            end
        end
        % differentiation
        varargout=diff(varargin)
        % unary functions
        varargout=sign(u)
        varargout=erf(u)
        varargout=uminus(u)
        varargout=sqrt(u)
        varargout=abs(u)
        varargout=log(u)
        varargout=log10(u)
        varargout=exp(u)
        varargout=cos(u)
        varargout=sin(u)
        varargout=tan(u)
        varargout=acos(u)
        varargout=asin(u)
        varargout=atan(u)
        varargout=cosh(u)
        varargout=sinh(u)
        varargout=tanh(u)
        varargout=acosh(u)
        varargout=asinh(u)
        varargout=atanh(u)
        % binary functions
        varargout=plus(u,v)
        varargout=minus(u,v)
        varargout=mtimes(u,v)
        varargout=mrdivide(u,v)
        varargout=mpower(u,v)
        varargout=max(u,v)
        varargout=min(u,v)
        varargout=gt(u,v)
        varargout=ge(u,v)
        varargout=eq(u,v)
        varargout=lt(u,v)
        varargout=le(u,v)
        varargout=ne(u,v)
        varargout=kron(u,v)
        % trinary functions
        varargout=normalpdf(u,v,w)
        varargout=normalcdf(u,v,w)
        varargout=if_then_else(u,v,w)
        varargout=if_elseif(varargin)
        % utility functions
        varargout=load_varlist(varargin)
        varargout=is_zero(varargin)
        varargout=is_one(varargin)
        varargout=isnumeric(varargin)
        varargout=recompose(varargin)
        %         varargout=char(varargin)
    end
    methods(Access=private)
        varargout=commit(varargin)
        varargout=compose_derivatives(varargin)
    end
    methods(Static)
        varargout=get_incidence(varargin)
        varargout=multinary_operation(varargin)
        varargout=commute(varargin)
        varargout=differentiate(varargin)
        varargout=push(varargin)
        varargout=arguments(varargin)
        varargout=system(varargin)
        varargout=swap_references(varargin)
%%%%%%%%%%%        varargout=transition_matrices2transition_matrix(varargin);
        varargout=extend_differentiation_list(varargin);
    end
    methods(Static,Access=private)
        varargout=initialize_differentiation_session(varargin)
        varargout=close_differentiation_session(varargin)
        varargout=equation2rise_sym(varargin)
        varargout=print(varargin)
        varargout=trim(varargin)
    end
end