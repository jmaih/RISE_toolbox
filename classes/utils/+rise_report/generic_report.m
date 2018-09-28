classdef generic_report < handle
    % Internal Object -- Not intended to be used directly by users
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
        
        Caption
        
        id
        
    end
    
    properties(Dependent)
        
        title
        
        subtitle
        
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
        
        function tt=get.title(obj)
            
            tt='';
            
            if ischar(obj.Caption)
                
                tt=obj.Caption;
                
            elseif iscell(obj.Caption)
                
                tt=obj.Caption{1};
                
            end
            
        end
        
        function stt=get.subtitle(obj)
            
            stt='';
            
            if iscell(obj.Caption)
                
                try
                    
                    stt=obj.Caption{2};
                    
                catch
                    
                    stt='';
                    
                end
                
            end
            
        end
        
        function C = reformat_text(obj,usertxt)
                        
            br = sprintf('\n');
            
            C = '';
            
            titre=obj.title;
                            
            if ~isempty(titre)
                
                C = [C,'\begin{tabular}{c}', br ];
                
                C = [C,format_title(), br ];
                
                C = [C,'\end{tabular}', br ];
                
            end
            
            if isempty(usertxt)
                
                return
                
            end
            
            if verbatim()
                
                C = [C,'\begin{verbatim}'];
                
            elseif ~centering()
                
                C = [C,'\begin{flushleft}'];
                
            end
            
            C = [C, br, usertxt];
            
            if verbatim()
                
                C = [C, br, '\end{verbatim}'];
                
            elseif ~centering()
                
                C = [C, br, '\end{flushleft}'];
                
            end
            
            C = [C,footnotetext()];
            
            function C = footnotetext()
                
                C = '';
                
                if any(strcmp(properties(obj),'footnote'))%isfield(obj,'footnote')
                    
                    footnote=obj.footnote;
                    
                    if ischar(footnote)
                        
                        footnote=cellstr(footnote);
                        
                    end
                    
                    if ~isempty(footnote)
                        
                        C = [footnote{:}];
                        
                    end
                    
                end
                
            end
            
            function C=format_title()
                
                C=rise_report.generic_report.reprocess(titre);
                
                if any(strcmp(properties(obj),'typeface'))%<---isfield(obj,'typeface')
                    
                    typeface=obj.typeface;
                    
                    if ischar(typeface)
                        
                        typeface=cellstr(typeface);
                        
                    end
                    
                    for it=1:numel(typeface)
                        
                        C=[typeface{it},'{',C,'}']; %#ok<AGROW>
                        
                    end
                    
                end
                
            end
            
            function flag=verbatim()
                
                flag=isfield(obj,'verbatim') && obj.verbatim;
                
            end
            
            function flag=centering()
                
                flag=isfield(obj,'centering') && obj.centering;
                
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
            
            if iscell(xin)
                
                for ii=1:numel(xin)
                    
                    xin{ii}=rise_report.generic_report.reprocess(xin{ii},precision);
                    
                end
                
                return
                
            end
            
            if isnumeric(xin)
            
                xin=num2str(xin,precision);
            
            elseif ischar(xin)
            
                if any(xin=='$')% ||sum(xin=='_')<=1 %&& ~isempty(strfind(xin,'ensuremath'))
                
                    return
                
                end
                
                repList={'_','%','&'};
                
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