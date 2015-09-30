classdef coef
    % coef class for setting linear restrictions on model parameters
    %
    % methods
    % --------
    %
    % - [block_exogenous](coef/block_exogenous)
    % - [coef](coef/coef)
    % - [linear_restrictions](coef/linear_restrictions)
    % - [minus](coef/minus)
    % - [mrdivide](coef/mrdivide)
    % - [mtimes](coef/mtimes)
    % - [plus](coef/plus)
    %
    % properties
    % -----------
    %
    % - [par_name] -
    properties
        % place holder for the parameter name, which is either
        %    - a string: this assumes that the name of the parameter is
        %    fully known in advance
        %    - a 3-element cell {eqtn,var_name,lag}: this information will
        %    be used later on to construct the string name of the
        %    corresponding parameter. The first element is the equation
        %    number, the second is the name of the endogenous variable and
        %    the third is the lag position.
        %    -  a 5-element cell {eqtn,var_name,lag,chain_name,state}.  Two
        %    elements are added to the previous case and they are needed in
        %    case the parameter is switching. So, the fourth element
        %    corresponds to the name of the markov chain controling the
        %    parameter and the fifth element is the state of the chain.
        par_name
    end
    properties(Access=protected)%,Hidden=true
        % place holder for function, which can only be one of the
        % following: plus, minus, mtimes, mrdivide
        func
        % place holder for the arguments of the function
        args
    end
    methods
        varargout=plus(varargin)
        varargout=minus(varargin)
        varargout=mtimes(varargin)
        varargout=mrdivide(varargin)
        varargout=linear_restrictions(varargin)
        function obj=coef(varargin)
            % there are 3 syntaxes for constructing a coef object
            % 1-) obj=coef('name')
            % 2-) obj=coef(eqtn,vname,lag)
            % 3-) obj=coef(eqtn,vname,lag,chain_name,state)
            n=nargin;
            if n
                if isa(varargin{1},mfilename)
                    obj=varargin{1};
                else
                    if n==1
                        if ~ischar(varargin{1})
                            error('with one argument, coef expects a string or a coef')
                        end
                        obj.par_name=varargin{1};
                    else
                        obj.par_name=varargin;
                    end
                end
            end
        end
    end
    methods(Static)
        varargout=block_exogenous(varargin)
    end
end