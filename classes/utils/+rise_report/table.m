classdef table < rise_report.generic_report
    properties
        title=''
        longtable=false
        % cell elements to be written in a table
        log={};
        % number of decimals
        precision=4
        numbering=true
    end
    properties(Dependent)
        batch
    end
    properties(Hidden)
        table_number
    end
    methods
        function obj=table(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b={};
            if ~isempty(obj.log)
                if ~isempty(obj.title)
                    best_title=@(x)rise_report.generic_report.best_title(x);
                    titel=obj.title;
                    if obj.numbering
                        titel=sprintf('Table \\# %0.0f: %s',obj.table_number,titel);
                    end
                    b=[b
                        best_title(titel)];
                end
                [nrows,ncols]=size(obj.log);
                AllBatch=cell(nrows,1);
                reprocess=@(x)rise_report.generic_report.reprocess(x,obj.precision);
                for irow=1:nrows
                    str='';
                    for icol=1:ncols
                        tmp=reprocess(obj.log{irow,icol});
                        % the leftmost elements are probably strings,
                        % aligned them left 
                        if icol==1
                            tmp=['\multicolumn{1}{l}{',tmp,'}']; %#ok<*AGROW>
                        end
                        str=[str,tmp];
                        if icol<ncols 
                            % skip to next column
                            str=[str,' &'];
                        end
                    end 
                    % end of line
                    str=[str,' \\'];
                    AllBatch{irow}=str;
                end
                theHeader=AllBatch{1};
                AllBatch=AllBatch(2:end);
                tablestyle='table';
                if obj.longtable
                    tablestyle=['long',tablestyle];
                    b=[b;{
                        ['\begin{',tablestyle,'}[H]{',repmat('r',1,ncols),'}']
                        }];
                else
                    b=[b;{
                        ['\begin{',tablestyle,'}[H] \centering']
                        ['\begin{tabular}{',repmat('r',1,ncols),'}']
                        }];
                end
                b=[b;{
                    '\hline\hline'
                    theHeader
                    '\hline'
                    }];
                if obj.longtable
                    b=[b;{
                        '\endfirsthead'
                        ['\multicolumn{',int2str(ncols),'}{c}{ -- \textit{Continued from previous page}} \\']
                        '\hline\hline'
                        theHeader
                        '\hline'
                        '\endhead'
                        ['\hline \multicolumn{',int2str(ncols),'}{r}{\textit{Continued on next page}} \\']
                        '\endfoot'
                        '\hline\hline'
                        '\endlastfoot'
                        }];
                end
                b=[b
                    AllBatch];
                if ~obj.longtable
                    b=[b
                        {'\hline\hline'
                        '\end{tabular}'}
                        ];
                end
                b=[b;{['\end{',tablestyle,'}']}];
            end
        end
    end
end