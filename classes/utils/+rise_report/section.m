classdef section < rise_report.generic_report
    % section reporting section
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.section/addlistener)
    % - [best_title](rise_report.section/best_title)
    % - [delete](rise_report.section/delete)
    % - [eq](rise_report.section/eq)
    % - [findobj](rise_report.section/findobj)
    % - [findprop](rise_report.section/findprop)
    % - [ge](rise_report.section/ge)
    % - [gt](rise_report.section/gt)
    % - [isvalid](rise_report.section/isvalid)
    % - [le](rise_report.section/le)
    % - [lt](rise_report.section/lt)
    % - [ne](rise_report.section/ne)
    % - [notify](rise_report.section/notify)
    % - [reprocess](rise_report.section/reprocess)
    % - [section](rise_report.section/section)
    % - [write](rise_report.section/write)
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
        function obj=section(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title,obj.numbering);
        end
    end
end
