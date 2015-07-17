classdef pagebreak < rise_report.generic_report
    % pagebreak report page break object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.pagebreak/addlistener)
    % - [best_title](rise_report.pagebreak/best_title)
    % - [delete](rise_report.pagebreak/delete)
    % - [eq](rise_report.pagebreak/eq)
    % - [findobj](rise_report.pagebreak/findobj)
    % - [findprop](rise_report.pagebreak/findprop)
    % - [ge](rise_report.pagebreak/ge)
    % - [gt](rise_report.pagebreak/gt)
    % - [isvalid](rise_report.pagebreak/isvalid)
    % - [le](rise_report.pagebreak/le)
    % - [lt](rise_report.pagebreak/lt)
    % - [ne](rise_report.pagebreak/ne)
    % - [notify](rise_report.pagebreak/notify)
    % - [pagebreak](rise_report.pagebreak/pagebreak)
    % - [reprocess](rise_report.pagebreak/reprocess)
    % - [write](rise_report.pagebreak/write)
    %
    % properties
    % -----------
    %
    % - [batch] -
    % - [id] -
    properties
        batch
    end
    methods
        function obj=pagebreak()
            obj.batch={' ';'\pagebreak'};
        end
    end
end
