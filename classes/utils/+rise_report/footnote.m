classdef footnote < rise_report.generic_report
    % footnote report footnote object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.footnote/addlistener)
    % - [best_title](rise_report.footnote/best_title)
    % - [delete](rise_report.footnote/delete)
    % - [eq](rise_report.footnote/eq)
    % - [findobj](rise_report.footnote/findobj)
    % - [findprop](rise_report.footnote/findprop)
    % - [footnote](rise_report.footnote/footnote)
    % - [ge](rise_report.footnote/ge)
    % - [gt](rise_report.footnote/gt)
    % - [isvalid](rise_report.footnote/isvalid)
    % - [le](rise_report.footnote/le)
    % - [lt](rise_report.footnote/lt)
    % - [ne](rise_report.footnote/ne)
    % - [notify](rise_report.footnote/notify)
    % - [reprocess](rise_report.footnote/reprocess)
    % - [write](rise_report.footnote/write)
    %
    % properties
    % -----------
    %
    % - [text] -
    % - [batch] -
    % - [id] -
    properties
        text=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=footnote(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            if ischar(obj.text)
                obj.text=cellstr(obj.text);
            end
        end
        function b = get.batch(obj)
            b=['\footnote{'
                obj.text(:)
                '}'];
        end
    end
end