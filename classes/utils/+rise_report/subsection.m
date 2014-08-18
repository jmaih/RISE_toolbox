classdef subsection < rise_report.generic_report
    properties
        title=''
        numbering=true
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=subsection(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title,obj.numbering);
        end
    end
end
