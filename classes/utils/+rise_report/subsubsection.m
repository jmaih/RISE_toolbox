classdef subsubsection < rise_report.generic_report
    % subsubsection report subsubsection object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.subsubsection/addlistener)
    % - [best_title](rise_report.subsubsection/best_title)
    % - [delete](rise_report.subsubsection/delete)
    % - [eq](rise_report.subsubsection/eq)
    % - [findobj](rise_report.subsubsection/findobj)
    % - [findprop](rise_report.subsubsection/findprop)
    % - [ge](rise_report.subsubsection/ge)
    % - [gt](rise_report.subsubsection/gt)
    % - [isvalid](rise_report.subsubsection/isvalid)
    % - [le](rise_report.subsubsection/le)
    % - [lt](rise_report.subsubsection/lt)
    % - [ne](rise_report.subsubsection/ne)
    % - [notify](rise_report.subsubsection/notify)
    % - [reprocess](rise_report.subsubsection/reprocess)
    % - [subsubsection](rise_report.subsubsection/subsubsection)
    % - [write](rise_report.subsubsection/write)
    %
    % properties
    % -----------
    %
    % - [title] -
    % - [numbering] -
    % - [batch] -
    % - [id] -
    properties
        title=''
        numbering=true
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=subsubsection(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title,obj.numbering);
        end
    end
end
