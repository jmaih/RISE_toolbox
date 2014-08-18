classdef footnote < rise_report.generic_reportdb
    properties
        text=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=footnote(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            if ischar(obj.text)
                obj.text=cellstr(obj.text);
            end
        end
        function b = get.batch(obj)
            b=['\footnote{'
                obj.text(:)
                '}'];
        end
    end
end