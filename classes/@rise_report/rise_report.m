classdef rise_report < handle
    properties
        documentclass='article'
        orientation='landscape'
        points='12pt'
        papersize='letterpaper'
        packages={'graphicx','amsmath','geometry','amsfonts',...
            'color','hyperref','longtable','float'}
        report_name='RiseReport'
        graphicspath=''
        pagestyle='myheadings'
        titlepage=struct('title','','date','','address','','author','',...
            'email','','abstract','')
    end
    properties(Hidden=true)
        compiler=getappdata(0,'rise_pdflatex')
    end
    properties(SetAccess = private,Hidden=true)
        script=cell(1000,1)
        ncells=1000;
        line_count=0
        lower_bound=100
        iter=1000
        tmpdir_flag=false;
        tmpdir='';
    end
    methods(Access=private)
        function record(obj,item)
            if ischar(item)
                item=cellstr(item);
            end
            if ~iscellstr(item)
                error('input must be a cellstr')
            end
            item=item(:);
            nrows=size(item,1);
            while obj.ncells-nrows-obj.line_count<obj.lower_bound
                obj.script(obj.line_count+(1:obj.iter))={};
                obj.ncells=obj.ncells+obj.iter;
            end
            obj.script(obj.line_count+(1:nrows))=item;
            % add some space between records for lisibility and for
            % debugging purposes
            obj.script(obj.line_count+nrows+1)={' '}; 
            obj.line_count=obj.line_count+nrows+1;
        end
    end
    methods
        function obj=rise_report(varargin)
            nargs=length(varargin);
            if nargs==0
                return
            end
            if isempty(obj.compiler)
                error([mfilename,':: cannot generate a report, MIKTEX was not found earlier'])
            end
            dummy=rise_report();
            default_field_map=fieldnames(dummy);
            default_field_map=default_field_map(:);
            nn=numel(default_field_map);
            default_report=struct();
            for ii=1:nn
                default_report.(default_field_map{ii})=obj.(default_field_map{ii});
            end
            is_update=isa(varargin{1},'rise_report');
            if is_update
                obj=varargin{1};
                if isempty(obj)
                    for irow=1:size(default_field_map,1)
                        disp(default_report.(default_field_map{irow,1}))
                    end
                    return
                end
                nargs=nargs-1;
                varargin=varargin(2:end);
                default_field_map=[default_field_map,num2cell(false(nn,1))];
            else
                default_field_map=[default_field_map,num2cell(true(nn,1))];
            end
            if rem(nargs,2)
                error('arguments must come in pairs')
            end
            for ii=1:2:nargs
                field=varargin{ii};
                loc=find(strcmp(field,default_field_map(:,1)));
                if isempty(loc)
                    error([field,':: unrecognized as a valid field for rise_report'])
                end
                value=varargin{ii+1};
                switch field
                    case 'documentclass'
                        default_field_map{loc,2}=true;
                        if iscellstr(value)
                            value=char(value);
                        end
                        if ~ischar(value)||~ismember(value,{'article','report'})
                            error('argument of documentclass must be a ''article'' or ''report'' ')
                        end
                        default_report.(field)=value;
                    case 'orientation'
                        default_field_map{loc,2}=true;
                        if ischar(value)
                            value=cellstr(value);
                        end
                        if ~iscellstr(value)||~ismember(value{1},{'portrait','landscape'})
                            error('argument of options must be a ''landscape'',''portrait''')
                        end
                        default_report.(field)=value{1};
                    case 'points'
                        default_field_map{loc,2}=true;
                        if ischar(value)
                            value=cellstr(value);
                        end
                        if ~iscellstr(value)||~ismember(value{1},{'10pt','11pt','12pt'})
                            error('argument of points must be a ''10pt'',''11pt'' or ''12pt''')
                        end
                        default_report.(field)=value{1};
                    case {'report_name','graphicspath'}
                        default_field_map{loc,2}=true;
                        if ~ischar(value)
                            error([field,' must be a char'])
                        end
                        default_report.(field)=value;
                    case 'pagestyle'
                        default_field_map{loc,2}=true;
                        if ~ismember(value,{'myheadings','empty','headings'})
                            error('argument of pagestyle must be a ''myheadings'',''empty'' or ''headings''')
                        end
                        default_report.(field)=value;
                    case 'papersize'
                        default_field_map{loc,2}=true;
                        if ~ismember(value,{'letterpaper','a4paper'})
                            error('argument of papersize must be a ''letterpaper'' or ''a4paper''')
                        end
                        default_report.(field)=value;
                    case 'packages'
                        default_field_map{loc,2}=true;
                        if ischar(value)
                            value=cellstr(value);
                        end
                        for ival=1:numel(value)
                            v=value{ival};
                            remove=v(1)=='-';
                            if remove
                                v=v(2:end);
                                loc=strcmp(v,default_report.packages);
                                default_report.packages(loc)=[];
                            else
                                default_report.packages=union(default_report.packages,v);
                            end
                        end
                    case 'titlepage'
                        default_field_map{loc,2}=true;
                        if ~isstruct(value)
                            error('the argument to titlepage must be a struct')
                        end
                        fields=fieldnames(value);
                        valid_fields=fieldnames(default_report.titlepage);
                        for ival=1:numel(fields)
                            ff=fields{ival};
                            if ~ismember(ff,valid_fields)
                                error([ff,' is not a valid field of titlepage'])
                            end
                            switch ff
                                case 'abstract'
                                    if ~isempty(value.(ff))
                                        if ischar(value.(ff))
                                            value.(ff)=cellstr(value.(ff));
                                        end
                                        default_report.titlepage.(ff)=value.(ff);
                                    end
                                case 'title'
                                    default_report.titlepage.(ff)=value.(ff);
                                case {'author','address','email'}
                                    if iscellstr(value.(ff))
                                        tmp=value.(ff){1};
                                        for jj=2:numel(value.(ff))
                                            tmp=[tmp,' \and ',value.(ff){jj}];
                                        end
                                        value.(ff)=tmp;
                                    end
                                    default_report.titlepage.(ff)=...
                                        regexp(value.(ff),'\\and(?!\w+)','split');
                                case 'date'
                                    default_report.titlepage.(ff)=value.(ff);
                                otherwise
                                    error('field not propertly handled in titlepage. Please contact junior.maih@gmail.com')
                            end
                        end
                    otherwise
                        error([field,':: is not handled. please report this problem to junior.maih@gmail.com'])
                end
            end
            for irow=1:size(default_field_map,1)
                if default_field_map{irow,2}
                    obj.(default_field_map{irow,1})=default_report.(default_field_map{irow,1});
                end
            end
        end
        %         %-----------------------------
        %         function matrix(obj,varargin)
        %             record(obj,varargin)
        %         end
        %         function equation(obj,varargin)
        %             %\begin{quotation}
        %             %\end{quotation}
        %             record(obj,varargin)
        %         end
        function subsubsection(obj,section_name,text)
            mysection={['\subsubsection{',char(section_name),'}']};
            if nargin>2
                if ischar(text)
                    text=cellstr(text);
                end
                mysection=[mysection;text(:)];
            end
            record(obj,mysection)
        end
        function subsection(obj,section_name,text)
            mysection={['\subsection{',char(section_name),'}']};
            if nargin>2
                if ischar(text)
                    text=cellstr(text);
                end
                mysection=[mysection;text(:)];
            end
            record(obj,mysection)
        end
        function section(obj,section_name,text)
            mysection={['\section{',char(section_name),'}']};
            if nargin>2
                if ischar(text)
                    text=cellstr(text);
                end
                mysection=[mysection;text(:)];
            end
            record(obj,mysection)
        end
        %         %-----------------------------
        function thisobj=table(obj,table_struct)
            default_table=struct('longtable',false,'title','',...
                'table',[],'caption','');
            if isempty(obj)
                thisobj=default_table;
                return
            end
            new_table=mysetfield(default_table,table_struct);
            myscript={};
            if ~isempty(new_table.title)
                myscript=[myscript;best_title(new_table.title)];
            end
            [nrows,ncols]=size(new_table.table);
            AllBatch=cell(nrows,1);
            for irow=1:nrows
                str='';
                for icol=1:ncols
                    tmp=reprocess(new_table.table{irow,icol});
                    if icol==1
                        tmp=['\multicolumn{1}{l}{',tmp,'}'];
                    end
                    str=[str,tmp];
                    if icol<ncols
                        str=[str,' &'];
                    end
                end
                str=[str,' \\'];
                AllBatch{irow}=str;
            end
            
            theHeader=AllBatch{1};
            AllBatch=AllBatch(2:end);
            tablestyle='table';
            if new_table.longtable
                tablestyle=['long',tablestyle];
                myscript=[myscript;{
                    ['\begin{',tablestyle,'}[H]{',repmat('r',1,ncols),'}']
                    '\hline\hline'
                    theHeader
                    '\hline'
                    }];
            else
                myscript=[myscript;{
                    ['\begin{',tablestyle,'}[H] \centering']
                    ['\begin{tabular}{',repmat('r',1,ncols),'}']
                    '\hline\hline'
                    theHeader
                    '\hline'
                    }];
            end
            if new_table.longtable
                myscript=[myscript;{
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
            myscript=[myscript
                AllBatch];
            if ~new_table.longtable
                myscript=[myscript
                    {'\hline\hline'
                    '\end{tabular}'}
                    ];
            end
            myscript=[myscript;{['\end{',tablestyle,'}']}];
            record(obj,myscript)
            thisobj=obj;
        end
        function thisobj=figure(obj,figure_struct)
            default_figure=struct('name','',...
                'title','',...
                'caption','','angle',0,'scale',0.85);
            if isempty(obj)
                thisobj=default_figure;
                return
            end
            new_figure=mysetfield(default_figure,figure_struct); %
            if ~isempty(new_figure.name)
                tmpname=new_figure.name;
                angle=new_figure.angle;
                scale=new_figure.scale;
                if ishandle(tmpname)
                    if ~obj.tmpdir_flag
                        obj.tmpdir=tempname(pwd);
                        mkdir(obj.tmpdir);
                        obj.tmpdir_flag=true;
                    end
                    new_handle=tmpname;
                    tmpname=tempname(obj.tmpdir);
                    angle=rise_saveaspdf(new_handle,tmpname);
                end
                tmpname=strrep(tmpname,'\','/');
                tmpname=strrep(tmpname,'.pdf','');%,'.pdf'
                % quotes and spaces may still not work if the extension .pdf is
                % added
                if any(isspace(tmpname))
                    tmpname=['"',tmpname,'"'];
                end
                myfigname=reprocess(new_figure.title);
%                 newscript={ %'\newpage'
%                     '\begin{tabular}[H]{@{\hspace*{-3pt}}c@{}}'
%                     ['\multicolumn{1}{c}{\large\bfseries ',myfigname,'}\\']
%                     ['\raisebox{10pt}{\includegraphics[scale=',num2str(scale),',angle=',...
%                     num2str(angle),']{',tmpname,'}}']
%                     '\end{tabular}'
%                     };
                newscript={ %'\newpage'
                    '\begin{tabular}[H]{@{\hspace*{-3pt}}c@{}}'
                    ['\multicolumn{1}{c}{\large\bfseries ',myfigname,'}\\']
                    ['\centerline{\includegraphics[scale=',num2str(scale),',angle=',...
                    num2str(angle),']{',tmpname,'}}']
                    '\end{tabular}'
                    };
                record(obj,newscript)
                thisobj=obj;
            end
        end
        function publish(obj,varargin)
            nargs=length(varargin);
            if nargs
                obj=rise_report(obj,varargin{:});
                if ~isempty(obj)
                    publish(obj)
                end
                return
            end
            if exist([obj.report_name,'.pdf'],'file')
                delete([obj.report_name,'.pdf'])
            end
            fclose('all');
            fid=fopen([obj.report_name,'.tex'],'w');
            add_preamble();
            add_title_page();
            for iline=1:obj.line_count
                fprintf(fid,'%s \n',obj.script{iline});
            end
            add_finishing();
            fclose(fid);
            
            thisString=[obj.compiler,' ',obj.report_name];
            retcode=system(thisString);
            if ~retcode
                % run again in order to make sure the references are shown
                system(thisString);
                system(thisString);
            end
            useless_extensions={'.log','.bbl','.blg','.aux','.*.bak'};
            for iext=1:numel(useless_extensions)
                file2delete=[obj.report_name,useless_extensions{iext}];
                if exist(file2delete,'file')
                    delete(file2delete)
                end
            end
            if obj.tmpdir_flag
                rmdir(obj.tmpdir,'s'); % use s to remove the contents as well otherwise the operation will not be successful
                obj.tmpdir_flag=false;
            end
            function add_finishing()
                fprintf(fid,'\n%s \n','\end{document}');
            end
            
            function add_preamble()
                preamble_instructions={
                    mydocumentclass()
                    myusepackage()%[centering]
                    '\numberwithin{equation}{section}'
                    '\begin{document}'
                    mypagestyle()
                    };
                if ~isempty(obj.graphicspath)
                    preamble_instructions=[preamble_instructions
                        ['\graphicspath{{',obj.graphicspath,'}}']];
                end
                for jj=1:numel(preamble_instructions)
                    fprintf(fid,'%s \n',preamble_instructions{jj});
                end
                fprintf(fid,'%s \n\n','% end of preamble');
                function out=mypagestyle()
                    out=['\pagestyle{',obj.pagestyle,'}'];
                end
                function out=myusepackage()
                    out='\usepackage{';
                    if ~isempty(obj.packages)
                        out=[out,obj.packages{1}];
                        for iopt=2:numel(obj.packages)
                            out=[out,',',obj.packages{iopt}]; %#ok<*AGROW>
                        end
                    end
                    out=[out,'}'];
                end
                function out=mydocumentclass()
                    out=['\documentclass[',obj.points,',',obj.orientation,']','{',obj.documentclass,'}'];
                end
            end
            
            function add_title_page()
                if isempty(obj.titlepage.title)
                    return
                end
                AfterAuthor=' \\';
                new_title=obj.titlepage;
                fprintf(fid,'%s \n',['\title{',reprocess(new_title.title),'}']);
                [authors,address,email]=collect_attributes(new_title,...
                    'author','address','email');
                flag_add=numel(address)==numel(authors);
                flag_email=numel(email)==numel(authors);
                fprintf(fid,'%s \n','\author{');
                for iaut=1:numel(authors)
                    thisAuthor=authors{iaut};
                    if flag_add
                        thisAuthor=[thisAuthor,AfterAuthor];
                        fprintf(fid,'%s \n',thisAuthor);
                        thisAddress=address{iaut};
                        if flag_email
                            thisAddress=[thisAddress,AfterAuthor];
                            fprintf(fid,'%s \n',thisAddress);
                            thisEmail=['\texttt{',email{iaut},'}'];
                            fprintf(fid,'%s \n',thisEmail);
                        else
                            fprintf(fid,'%s \n',thisAddress);
                        end
                    else
                        if flag_email
                            thisAuthor=[thisAuthor,AfterAuthor];
                            fprintf(fid,'%s \n',thisAuthor);
                            thisEmail=['\texttt{',email{iaut},'}'];
                            fprintf(fid,'%s \n',thisEmail);
                        else
                            fprintf(fid,'%s \n',thisAuthor);
                        end
                    end
                    if iaut<numel(authors)
                        fprintf(fid,'%s \n','\and');
                    end
                end
                fprintf(fid,'%s \n','}');
                if ~isempty(new_title.date)
                    fprintf(fid,'%s \n',['\date{',new_title.date,'}']);
                end
                fprintf(fid,'%s \n','\pagestyle{myheadings}');
                fprintf(fid,'%s \n',['\markright{',reprocess(new_title.title),'\hfill ',new_title.date,'\hfill}']);
                fprintf(fid,'%s \n','\thispagestyle{empty}');
                fprintf(fid,'%s \n','\maketitle');
                if ~isempty(new_title.abstract)
                    fprintf(fid,'%s \n','\begin{abstract}');
                    for irow=1:numel(new_title.abstract)
                        fprintf(fid,'%s \n',new_title.abstract{irow});
                    end
                    fprintf(fid,'%s \n','\end{abstract}');
                end
                function varargout=collect_attributes(new_title,varargin)
                    varargout=varargin;
                    for ii=1:length(varargin)
                        varargout{ii}=[];
                        if ~isempty(new_title.(varargin{ii}))
                            if ischar(new_title.(varargin{ii}))
                                new_title.(varargin{ii})=cellstr(new_title.(varargin{ii}));
                            end
                            varargout{ii}=new_title.(varargin{ii});
                        end
                    end
                end
            end
        end
        function report(obj,varargin)
            publish(obj,varargin)
        end
        function include(obj,filename)
            record(obj,['\input{',filename,'}']);
        end
        function verbatim(obj,thisstruct)
            default=struct('title','','list',{{}});
            if isempty(obj)
                disp(default)
                return
            end
            if ~isstruct(thisstruct)
                error('input must be a structure')
            end
            default=mysetfield(default,thisstruct);
             myscript={};
            if ~isempty(default.title)
                myscript=[myscript;best_title(default.title)];
            end
            list=default.list;
            if ~isempty(list)
                if ischar(list)
                    list=cellstr(list);
                end
                record(obj,[myscript;
                    {'\begin{verbatim}'};
                    list(:);
                    {'\end{verbatim}'}])
            end
        end
        function text(obj,thisstruct)
            default=struct('title','','list',{{}});
            if isempty(obj)
                disp(default)
                return
            end
            if ~isstruct(thisstruct)
                error('input must be a structure')
            end
            default=mysetfield(defaut,thisstruct);
             myscript={};
            if ~isempty(default.title)
                myscript=[myscript
                    {'\begin{tabular}{c}'}
                    {['\multicolumn{1}{c}{\large\bfseries ',reprocess(default.title),'}']}
                    {'\end{tabular}'}
                    ];
            end
            list=default.list;
            if ~isempty(list)
                if ischar(list)
                    list=cellstr(list);
                end
                record(obj,[myscript;
                    list(:)])
            end
        end
        function paragraph(obj,varargin)
            text(obj,...
                struct('title','','list',{[' ';varargin(:)]}))
        end
        function quote(obj,quote)
            if ischar(quote)
                quote=cellstr(quote);
            end
            record(obj,['\begin{quote}'
                quote(:)
                '\end{quote}'])
        end
        function quotation(obj,quote)
            if ischar(quote)
                quote=cellstr(quote);
            end
            record(obj,['\begin{quotation}'
                quote(:)
                '\end{quotation}'])
        end
        function enumerate(obj,thisstruct)
            default=struct('title','','list',{{}});
            if isempty(obj)
                disp(default)
                return
            end
            if ~isstruct(thisstruct)
                error('input must be a structure')
            end
            default=mysetfield(defaut,thisstruct);
             myscript={};
            if ~isempty(default.title)
                myscript=[myscript;best_title(default.title)];
            end
            list=default.list;
            if ~isempty(list)
                if ischar(list)
                    list=cellstr(list);
                end
                for ii=1:numel(list)
                    list{ii}=['\item ',list{ii}];
                end
                myscript=[myscript;
                    {'\begin{enumerate}'}
                    list(:)
                    {'\end{enumerate}'}
                    ];
                record(obj,myscript)
            end
        end
        function itemize(obj,thisstruct)
            default=struct('title','','list',{{}});
            if isempty(obj)
                disp(default)
                return
            end
            if ~isstruct(thisstruct)
                error('input must be a structure')
            end
            default=mysetfield(defaut,thisstruct);
             myscript={};
            if ~isempty(default.title)
                myscript=[myscript;best_title(default.title)];
            end
            list=default.list;
            if ~isempty(list)
                if ischar(list)
                    list=cellstr(list);
                end
                for ii=1:numel(list)
                    list{ii}=['\item ',list{ii}];
                end
                myscript=[myscript;
                    {'\begin{itemize}'}
                    list(:)
                    {'\end{itemize}'}
                    ];
                record(obj,myscript)
            end
        end
        function clearpage(obj)
            record(obj,{'\clearpage'})
        end
        function newpage(obj)
            record(obj,{'\newpage'})
        end
        function pagebreak(obj)
            record(obj,{'\pagebreak'})
        end
    end
end

function tt=best_title(tit)
tt={'\begin{tabular}{c}';
['\multicolumn{1}{c}{\large\bfseries ',reprocess(tit),'}']
'\end{tabular}'};
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