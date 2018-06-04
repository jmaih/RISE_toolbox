classdef table < rise_report.generic_report
% Create a table
%
% ::
%
%    document.table('option_name',option_value);
%
% Args:
%
%    varargin: arguments need to come in pairs
%
%       - **log** : [cell] contents of the table
%       - **title** : [char] title of the table
%       - **longtable** : [true|{false}] whether long table
%       - **numbering** : [{true}|false] whether to number the table
%       - **precision** : [numeric|{4}] number of decimal places
%       - **header** : [{''}|char|cell] : Header of the table if not already included in the log
%
% Example:
%
% ::
%
%    table={
%        '','col2','col3','col4'
%        'row2',rand,rand,rand
%        'row3',rand,rand,rand
%        'row4',rand,rand,rand
%        'row5',rand,rand,rand
%        'row6',rand,rand,rand
%        'row7',rand,rand,rand
%        'row8',rand,rand,rand
%        };
%    document.table('title','title of the table', 'log',table)
%



% methods
% --------
%
% - [addlistener](rise_report.table/addlistener)
% - [best_title](rise_report.table/best_title)
% - [delete](rise_report.table/delete)
% - [eq](rise_report.table/eq)
% - [findobj](rise_report.table/findobj)
% - [findprop](rise_report.table/findprop)
% - [ge](rise_report.table/ge)
% - [gt](rise_report.table/gt)
% - [isvalid](rise_report.table/isvalid)
% - [le](rise_report.table/le)
% - [lt](rise_report.table/lt)
% - [ne](rise_report.table/ne)
% - [notify](rise_report.table/notify)
% - [reprocess](rise_report.table/reprocess)
% - [table](rise_report.table/table)
% - [write](rise_report.table/write)
%
% properties
% -----------
%
% - [title] -
% - [longtable] -
% - [log] -
% - [precision] -
% - [numbering] -
% - [batch] -
% - [id] -
properties
title=''
longtable=false
% cell elements to be written in a table
log={};
% number of decimals
precision=4
numbering=true
header={};
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
                theHeader=obj.header;
                if isempty(theHeader)
                    theHeader=AllBatch(1);
                    AllBatch=AllBatch(2:end);
                end
                tablestyle='table';
                if obj.longtable
                    tablestyle=['long',tablestyle];
                    b=[b;
                        ['\begin{',tablestyle,'}[H]{',repmat('r',1,ncols),'}']
                        ];
                else
                    b=[b
                        ['\begin{',tablestyle,'}[H] \centering']
                        ['\begin{tabular}{',repmat('r',1,ncols),'}']
                        ];
                end
                if ~isempty(obj.title)
                    titel=obj.title;
                    if obj.numbering
                        titel=sprintf('Table \\# %0.0f: %s',obj.table_number,titel);
                    end
                    titel=reprocess(titel);
                    b=[
                        b
                        ['\multicolumn{',int2str(ncols),'}{c}{\large\bfseries ',titel,'}\\']
                        ];
                end
                b=[b
                    '\hline\hline'
                    theHeader
                    '\hline'
                    ];
                if obj.longtable
                    b=[b;
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
                        ];
                end
                b=[b
                    AllBatch];
                if ~obj.longtable
                    b=[b
                        '\hline\hline'
                        '\end{tabular}'
                        ];
                end
                b=[b
                    ['\end{',tablestyle,'}']
                    ];
            end
        end
    end
end