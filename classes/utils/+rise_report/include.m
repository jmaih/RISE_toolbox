classdef include < rise_report.generic_report
    % include report include object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.include/addlistener)
    % - [best_title](rise_report.include/best_title)
    % - [delete](rise_report.include/delete)
    % - [eq](rise_report.include/eq)
    % - [findobj](rise_report.include/findobj)
    % - [findprop](rise_report.include/findprop)
    % - [ge](rise_report.include/ge)
    % - [gt](rise_report.include/gt)
    % - [include](rise_report.include/include)
    % - [isvalid](rise_report.include/isvalid)
    % - [le](rise_report.include/le)
    % - [lt](rise_report.include/lt)
    % - [ne](rise_report.include/ne)
    % - [notify](rise_report.include/notify)
    % - [reprocess](rise_report.include/reprocess)
    % - [write](rise_report.include/write)
    %
    % properties
    % -----------
    %
    % - [filename] -
    % - [batch] -
    % - [id] -
    properties
        filename=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=include(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b={};
            if ~isempty(obj.filename)
                b={['\input{',obj.filename,'}']};
            end
        end
    end
end
