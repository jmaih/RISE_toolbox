classdef clearpage < rise_report.generic_report
    % clearpage report clear page object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.clearpage/addlistener)
    % - [best_title](rise_report.clearpage/best_title)
    % - [clearpage](rise_report.clearpage/clearpage)
    % - [delete](rise_report.clearpage/delete)
    % - [eq](rise_report.clearpage/eq)
    % - [findobj](rise_report.clearpage/findobj)
    % - [findprop](rise_report.clearpage/findprop)
    % - [ge](rise_report.clearpage/ge)
    % - [gt](rise_report.clearpage/gt)
    % - [isvalid](rise_report.clearpage/isvalid)
    % - [le](rise_report.clearpage/le)
    % - [lt](rise_report.clearpage/lt)
    % - [ne](rise_report.clearpage/ne)
    % - [notify](rise_report.clearpage/notify)
    % - [reprocess](rise_report.clearpage/reprocess)
    % - [write](rise_report.clearpage/write)
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
        function obj=clearpage()
            obj.batch={' ';'\clearpage'};
        end
    end
end
