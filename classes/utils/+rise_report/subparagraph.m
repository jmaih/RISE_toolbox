classdef subparagraph < rise_report.generic_report
    % subparagraph report subparagraph object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.subparagraph/addlistener)
    % - [best_title](rise_report.subparagraph/best_title)
    % - [delete](rise_report.subparagraph/delete)
    % - [eq](rise_report.subparagraph/eq)
    % - [findobj](rise_report.subparagraph/findobj)
    % - [findprop](rise_report.subparagraph/findprop)
    % - [ge](rise_report.subparagraph/ge)
    % - [gt](rise_report.subparagraph/gt)
    % - [isvalid](rise_report.subparagraph/isvalid)
    % - [le](rise_report.subparagraph/le)
    % - [lt](rise_report.subparagraph/lt)
    % - [ne](rise_report.subparagraph/ne)
    % - [notify](rise_report.subparagraph/notify)
    % - [reprocess](rise_report.subparagraph/reprocess)
    % - [subparagraph](rise_report.subparagraph/subparagraph)
    % - [write](rise_report.subparagraph/write)
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
        function obj=subparagraph(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title);
        end
    end
end
