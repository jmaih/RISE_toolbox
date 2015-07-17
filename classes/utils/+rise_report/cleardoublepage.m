classdef cleardoublepage < rise_report.generic_report
    % cleardoublepage report clear double page object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.cleardoublepage/addlistener)
    % - [best_title](rise_report.cleardoublepage/best_title)
    % - [cleardoublepage](rise_report.cleardoublepage/cleardoublepage)
    % - [delete](rise_report.cleardoublepage/delete)
    % - [eq](rise_report.cleardoublepage/eq)
    % - [findobj](rise_report.cleardoublepage/findobj)
    % - [findprop](rise_report.cleardoublepage/findprop)
    % - [ge](rise_report.cleardoublepage/ge)
    % - [gt](rise_report.cleardoublepage/gt)
    % - [isvalid](rise_report.cleardoublepage/isvalid)
    % - [le](rise_report.cleardoublepage/le)
    % - [lt](rise_report.cleardoublepage/lt)
    % - [ne](rise_report.cleardoublepage/ne)
    % - [notify](rise_report.cleardoublepage/notify)
    % - [reprocess](rise_report.cleardoublepage/reprocess)
    % - [write](rise_report.cleardoublepage/write)
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
        function obj=cleardoublepage()
            obj.batch={' ';'\cleardoublepage'};
        end
    end
end
