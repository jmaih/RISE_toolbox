classdef chapter < rise_report.generic_report
    % chapter report chapter object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.chapter/addlistener)
    % - [best_title](rise_report.chapter/best_title)
    % - [chapter](rise_report.chapter/chapter)
    % - [delete](rise_report.chapter/delete)
    % - [eq](rise_report.chapter/eq)
    % - [findobj](rise_report.chapter/findobj)
    % - [findprop](rise_report.chapter/findprop)
    % - [ge](rise_report.chapter/ge)
    % - [gt](rise_report.chapter/gt)
    % - [isvalid](rise_report.chapter/isvalid)
    % - [le](rise_report.chapter/le)
    % - [lt](rise_report.chapter/lt)
    % - [ne](rise_report.chapter/ne)
    % - [notify](rise_report.chapter/notify)
    % - [reprocess](rise_report.chapter/reprocess)
    % - [write](rise_report.chapter/write)
    %
    % properties
    % -----------
    %
    % - [title] -
    % - [tableOfContents] -
    % - [tableOfContentsTitle] -
    % - [numbering] -
    % - [batch] -
    % - [id] -
    properties
        title=''
        tableOfContents=false
        tableOfContentsTitle=''
        numbering=true % preface is a chapter and is not numbered
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=chapter(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title,obj.numbering,...
                obj.tableOfContentsTitle);
        end
    end
end
