classdef paragraph < rise_report.generic_report
    % paragraph report paragraph object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.paragraph/addlistener)
    % - [best_title](rise_report.paragraph/best_title)
    % - [delete](rise_report.paragraph/delete)
    % - [eq](rise_report.paragraph/eq)
    % - [findobj](rise_report.paragraph/findobj)
    % - [findprop](rise_report.paragraph/findprop)
    % - [ge](rise_report.paragraph/ge)
    % - [gt](rise_report.paragraph/gt)
    % - [isvalid](rise_report.paragraph/isvalid)
    % - [le](rise_report.paragraph/le)
    % - [lt](rise_report.paragraph/lt)
    % - [ne](rise_report.paragraph/ne)
    % - [notify](rise_report.paragraph/notify)
    % - [paragraph](rise_report.paragraph/paragraph)
    % - [reprocess](rise_report.paragraph/reprocess)
    % - [write](rise_report.paragraph/write)
    %
    % properties
    % -----------
    %
    % - [title] -
    % - [batch] -
    % - [id] -
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
