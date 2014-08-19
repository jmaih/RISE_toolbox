classdef clearpage < rise_report.generic_report
    properties
        batch
    end
    methods
        function obj=clearpage()
            obj.batch={' ';'\clearpage'};
        end
    end
end
