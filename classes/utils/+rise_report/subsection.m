classdef subsection < rise_report.generic_report
    % subection report section object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.subsection/addlistener)
    % - [best_title](rise_report.subsection/best_title)
    % - [delete](rise_report.subsection/delete)
    % - [eq](rise_report.subsection/eq)
    % - [findobj](rise_report.subsection/findobj)
    % - [findprop](rise_report.subsection/findprop)
    % - [ge](rise_report.subsection/ge)
    % - [gt](rise_report.subsection/gt)
    % - [isvalid](rise_report.subsection/isvalid)
    % - [le](rise_report.subsection/le)
    % - [lt](rise_report.subsection/lt)
    % - [ne](rise_report.subsection/ne)
    % - [notify](rise_report.subsection/notify)
    % - [reprocess](rise_report.subsection/reprocess)
    % - [subsection](rise_report.subsection/subsection)
    % - [write](rise_report.subsection/write)
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
        function obj=subsection(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title,obj.numbering);
        end
    end
end
