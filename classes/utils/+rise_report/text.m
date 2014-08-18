classdef text < rise_report.generic_report
    properties
        log={} % cellstring
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=text(log_)
            if ischar(log_)
                log_=cellstr(log_);
            end
            if ~iscellstr(log_)
                error([mfilename,':: input must be char or cellstr'])
            end
            obj.log=log_;
        end
        function b = get.batch(obj)
            b=batch_implementation(obj,obj.log(:));
        end 
    end
    methods(Access = private)
        function b=batch_implementation(~,b)
            % this function does nothing
        end
    end
end