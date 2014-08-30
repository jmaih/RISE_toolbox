classdef quotation < rise_report.generic_report
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