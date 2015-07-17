classdef verbatim < rise_report.text
    % verbatim report verbatim object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.verbatim/addlistener)
    % - [best_title](rise_report.verbatim/best_title)
    % - [delete](rise_report.verbatim/delete)
    % - [eq](rise_report.verbatim/eq)
    % - [findobj](rise_report.verbatim/findobj)
    % - [findprop](rise_report.verbatim/findprop)
    % - [ge](rise_report.verbatim/ge)
    % - [gt](rise_report.verbatim/gt)
    % - [isvalid](rise_report.verbatim/isvalid)
    % - [le](rise_report.verbatim/le)
    % - [lt](rise_report.verbatim/lt)
    % - [ne](rise_report.verbatim/ne)
    % - [notify](rise_report.verbatim/notify)
    % - [reprocess](rise_report.verbatim/reprocess)
    % - [verbatim](rise_report.verbatim/verbatim)
    % - [write](rise_report.verbatim/write)
    %
    % properties
    % -----------
    %
    % - [log] -
    % - [batch] -
    methods
        function obj=verbatim(varargin)
            obj=obj@rise_report.text(varargin{:});
        end
    end
    methods(Access = protected)
        function b=batch_implementation(~,b)
            b=['\begin{verbatim}';b(:);'\end{verbatim}'];
        end
    end
end
