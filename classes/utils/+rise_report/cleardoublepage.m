classdef cleardoublepage < rise_report.generic_report
    properties
        batch
    end
    methods
        function obj=cleardoublepage()
            obj.batch={'\cleardoublepage'};
        end
    end
end
