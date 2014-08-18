classdef pagebreak < rise_report.generic_report
    properties
        batch
    end
    methods
        function obj=pagebreak()
            obj.batch={'\pagebreak'};
        end
    end
end
