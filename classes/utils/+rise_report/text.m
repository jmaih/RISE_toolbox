classdef text < rise_report.generic_report
    % text report text object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.text/addlistener)
    % - [best_title](rise_report.text/best_title)
    % - [delete](rise_report.text/delete)
    % - [eq](rise_report.text/eq)
    % - [findobj](rise_report.text/findobj)
    % - [findprop](rise_report.text/findprop)
    % - [ge](rise_report.text/ge)
    % - [gt](rise_report.text/gt)
    % - [isvalid](rise_report.text/isvalid)
    % - [le](rise_report.text/le)
    % - [lt](rise_report.text/lt)
    % - [ne](rise_report.text/ne)
    % - [notify](rise_report.text/notify)
    % - [reprocess](rise_report.text/reprocess)
    % - [text](rise_report.text/text)
    % - [write](rise_report.text/write)
    %
    % properties
    % -----------
    %
    % - [log] -
    % - [batch] -
    % - [id] -
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
    methods(Access = protected)
        function b=batch_implementation(~,b)
            % this function does nothing
        end
    end
end