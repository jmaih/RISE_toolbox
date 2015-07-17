classdef quotation < rise_report.generic_report
    % quotation report quotation object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.quotation/addlistener)
    % - [best_title](rise_report.quotation/best_title)
    % - [delete](rise_report.quotation/delete)
    % - [eq](rise_report.quotation/eq)
    % - [findobj](rise_report.quotation/findobj)
    % - [findprop](rise_report.quotation/findprop)
    % - [ge](rise_report.quotation/ge)
    % - [gt](rise_report.quotation/gt)
    % - [isvalid](rise_report.quotation/isvalid)
    % - [le](rise_report.quotation/le)
    % - [lt](rise_report.quotation/lt)
    % - [ne](rise_report.quotation/ne)
    % - [notify](rise_report.quotation/notify)
    % - [quotation](rise_report.quotation/quotation)
    % - [reprocess](rise_report.quotation/reprocess)
    % - [write](rise_report.quotation/write)
    %
    % properties
    % -----------
    %
    % - [log] -
    % - [reference] -
    % - [batch] -
    % - [id] -
    properties
        log={} % cellstring
        reference
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=quotation(log_,ref)
            if nargin
                if ischar(log_)
                    log_=cellstr(log_);
                end
                if ~iscellstr(log_)
                    error([mfilename,':: input must be char or cellstr'])
                end
                obj.log=log_;
                if nargin>1
                    obj.reference=char(ref);
                end
            end
        end
        function b = get.batch(obj)
            b={};
            if ~isempty(obj.log)
                b=[b % <-- this ensures cell
                    %                     '\begin{tabular}[H]{@{\hspace*{-3pt}}c@{}}'
                    '\begin{quotation}'
                    obj.log(:)
                    '\end{quotation}'
                    ];
                if ~isempty(obj.reference)
                    b=[b
                        '\bigskip'
                        %                         ['\multicolumn{1}{c}{\textit\bfseries ',char(obj.reference),'}\\']
                        obj.reference(:)
                        ];
                end
                %                 b=[b
                %                     '\end{tabular}'
                %                     ];
            end
        end
    end
end