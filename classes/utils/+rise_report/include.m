classdef include < rise_report.generic_report
    properties
        filename=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=include(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b={};
            if ~isempty(obj.filename)
                b={['\input{',obj.filename,'}']};
            end
        end
    end
end
