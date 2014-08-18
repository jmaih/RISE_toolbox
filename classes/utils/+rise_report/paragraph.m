classdef paragraph < rise_report.generic_report
    properties
        title=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=paragraph(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title);
        end
    end
end
