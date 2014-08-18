classdef generic_report < handle
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
                xin=strrep(xin,'\_','LouisPergaud');
                xin=strrep(xin,'_','\_');
                xin=strrep(xin,'LouisPergaud','\_');
            end
            xin=strtrim(xin);
        end
    end
end