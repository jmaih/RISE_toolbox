classdef rise_equation
    properties(SetAccess=protected)
        dynamic
        static
        shadow_dynamic
        shadow_static
        shadow_balanced_growth_path
        tex_dynamic
        tex_static
        id
    end
    methods
        % set property utility
        function obj=set_properties(obj,varargin)
            nargs=length(varargin);
            if rem(nargs,2)~=0
                error([mfilename,':: arguments should enter by pairs'])
            end
            object_properties=properties(obj);
            for ii=nargs/2:-1:1
                propname=varargin{2*ii-1};
                if isempty(strcmp(propname,object_properties))
                    error([mfilename,...
                        ':: propname is not a valid property of class ',class(obj)])
                end
                propval=varargin{2*ii};
                obj.(propname)=propval;
            end
        end
        % constructor
        function obj=rise_equation(varargin)
            if nargin
                obj=set_properties(obj,varargin{:});
            end
        end
    end
end

