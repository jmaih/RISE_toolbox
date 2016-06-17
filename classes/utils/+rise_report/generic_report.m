classdef generic_report < handle
    % generic_report generic class for report objects
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.generic_report/addlistener)
    % - [best_title](rise_report.generic_report/best_title)
    % - [delete](rise_report.generic_report/delete)
    % - [eq](rise_report.generic_report/eq)
    % - [findobj](rise_report.generic_report/findobj)
    % - [findprop](rise_report.generic_report/findprop)
    % - [ge](rise_report.generic_report/ge)
    % - [generic_report](rise_report.generic_report/generic_report)
    % - [gt](rise_report.generic_report/gt)
    % - [isvalid](rise_report.generic_report/isvalid)
    % - [le](rise_report.generic_report/le)
    % - [lt](rise_report.generic_report/lt)
    % - [ne](rise_report.generic_report/ne)
    % - [notify](rise_report.generic_report/notify)
    % - [reprocess](rise_report.generic_report/reprocess)
    % - [write](rise_report.generic_report/write)
    %
    % properties
    % -----------
    %
    % - [id] -
    properties
        id
    end
    methods
        function write(obj,fid)
            b=obj.batch;
            environment='%-------------------------------------------------%';
            fprintf(fid,'%s\n',environment);
            fprintf(fid,'%s item # %0.0f\n','%',obj.id);
            fprintf(fid,'%s\n',environment);
            for ii=1:numel(b)
                fprintf(fid,'%s \n',b{ii});
            end
        end
    end
    methods(Static)
        function tt=best_title(tit)
            tt={
                '\begin{tabular}[c]{c}';
                ['\multicolumn{1}{c}{\large\bfseries ',...
                rise_report.generic_report.reprocess(tit),'}']
                '\end{tabular}'
                };
        end
        function xin=reprocess(xin,precision)
            if nargin<2
                precision=4;
            end
            if isnumeric(xin)
                xin=num2str(xin,precision);
            elseif ischar(xin)
                if any(xin=='$')% ||sum(xin=='_')<=1 %&& ~isempty(strfind(xin,'ensuremath'))
                    return
                end
                repList={'_','%'};
                for ilist=1:numel(repList)
                    % protect the already escaped
                    xin=strrep(xin,['\',repList{ilist}],'LouisPergaud');
                    % escape the un-escaped
                    xin=strrep(xin,repList{ilist},['\',repList{ilist}]);
                    % undo the protection
                    xin=strrep(xin,'LouisPergaud',['\',repList{ilist}]);
                end
            end
            xin=strtrim(xin);
        end
    end
end