classdef rise_report < handle
    properties
        documentclass='article'
        orientation='landscape'
        points='12pt'
        packages={'graphicx','amsmath','geometry','amsfonts',...
            'color','hyperref','longtable'}
        report_name='RiseReport'
        graphicspath=''
        pagestyle='myheadings'
        titlepage=struct('title','','date','','address','','author','','email','')
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
            %             for irow=1:nrows
            %                 obj.script(obj.line_count+irow)=item(irow);
            %             end
            obj.line_count=obj.line_count+nrows;
        end
    end
    methods
        function obj=rise_report(varargin)
            nargs=length(varargin);
            if nargs==0
                return
            end
            rise_pdflatex=getappdata(0,'rise_pdflatex');
            if ~rise_pdflatex
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
                                case 'title'
                                    default_report.titlepage.(ff)=value.(ff);
                                case {'author','address','email'}
                                    if iscellstr(value.(ff))
                                        tmp=value.(ff);
                                        tmp=cell2mat(strcat(tmp(:)',' &'));
                                        tmp=tmp(1:end-1);
                                        value.(ff)=tmp;
                                    end
                                    default_report.titlepage.(ff)=value.(ff);
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
        %         function model_equations(obj,model)
        %         end
        %         function model_solution(obj,model)
        %         end
        %         %-----------------------------
        function model_equations(obj,model_objects)
            new_equations=struct('title','Model equations',...
                'equations',{{model_objects(1).equations.dynamic}});
            if ischar(new_equations.equations)
                new_equations.equations=cellstr(new_equations.equations);
            end
            newscript={
                ['\section{',reprocess(new_equations.title),'}']
                ' \begin{verbatim}'	%\color{lightgray}
                };
            for ieq=1:numel(new_equations.equations)
                myeq=['EQ',int2str(ieq),': ',new_equations.equations{ieq}];
                equality=strfind(myeq,'=');
                if isempty(equality)
                    myeq=strrep(myeq,';','=0;');
                end
                if length(myeq)>89
                    neweq='';
                    strfill='';
                    while ~isempty(myeq)
                        [token,remain] = strtok(myeq,'+-/*');
                        strfill=[strfill,token];
                        if ~isempty(remain)
                            strfill=[strfill,remain(1)]; %#ok<*AGROW>
                        end
                        myeq=remain(2:end);
                        if length(strfill)>86
                            if ~isempty(remain(2:end))
                                strfill=[strfill,'...'];
                            end
                            neweq=char(neweq,strfill);
                            if ~isempty(remain(2:end))
                                strfill='	';
                            else
                                strfill='';
                            end
                        end
                    end
                    if ~isempty(strfill)
                        neweq=char(neweq,strfill);
                    end
                    myeq=neweq(2:end,:);
                    clear neweq
                end
                for isubeq=1:size(myeq,1)
                    newscript=[newscript;{myeq(isubeq,:)}];
                end
                newscript=[newscript;{' '}];
            end
            newscript=[newscript;{'\end{verbatim}\color{black}'}];
            text(obj,newscript);
        end
        function model_estimation_results(obj,model_objects)
            ncases=numel(model_objects);
            type_name='tex_name';
            parnames= {model_objects(1).estimated_parameters.(type_name)};
            ordered_names=sort(parnames);
            PALL=cell(1,ncases);
            for ic=1:ncases
                newnames={model_objects(ic).estimated_parameters.(type_name)};
                mode=num2cell(vertcat(model_objects(ic).estimated_parameters.mode));
                mode_std=num2cell(vertcat(model_objects(ic).estimated_parameters.mode_std));
                prior_prob=num2cell(vertcat(model_objects(ic).estimated_parameters.interval_probability));
                plb=num2cell(vertcat(model_objects(ic).estimated_parameters.plb));
                pub=num2cell(vertcat(model_objects(ic).estimated_parameters.pub));
                prior_distrib={model_objects(ic).estimated_parameters.distribution};
                PALL{ic}=[newnames',prior_distrib',prior_prob,plb,pub,mode,mode_std];
                if ~isequal(parnames,newnames)
                    ordered_names=union(ordered_names,newnames);
                end
            end
            nparams=numel(ordered_names);
            
            mytable=cell(nparams+1,5+ncases);
            model_names={model_objects.filename};
            mytable(1,:)=[{'parameter','Prior distr','Prior prob','low','high'},model_names];
            for iparam=1:nparams
                name_in=false;
                for imod=1:ncases
                    if isempty(PALL{imod})
                        mytable{iparam+1,5+imod}='--';
                    else
                        if ~name_in
                            param_info=PALL{imod}(1,:);
                            pname=param_info{1};
                            name_in=true;
                            mytable(iparam+1,1:5)=param_info(1:5);
                            mytable{iparam+1,5+imod}=param_info{6};
                            PALL{imod}(1,:)=[];
                            % add the dollars
                            mytable{iparam+1,1}=['$',mytable{iparam+1,1},'$'];
                            mytable{iparam+1,1}=strrep(mytable{iparam+1,1},'_','\_');
                            continue
                        end
                        loc=find(strcmp(pname,PALL{imod}(:,1)));
                        if isempty(loc)
                            mytable{iparam+1,5+imod}='--';
                        else
                            mytable{iparam+1,5+imod}=PALL{imod}{loc,6};
                            PALL{imod}(loc,:)=[];
                        end
                    end
                end
            end
            if ncases==1
                mytable=[mytable,[{'mode\_std'};num2cell(vertcat(model_objects.estimated_parameters.mode_std))]];
            end
            table_struct=struct('longtable',true,'table',{mytable},...
                'caption','Prior and Posterior mode',...
                'title','Estimation Results');
            table(obj,table_struct)
            % add estimation statistics
            estimation_statistics()
            function estimation_statistics()
                thisTable={
                    ''
                    'log-post:'
                    'log-lik:'
                    'log-prior:'
                    'log-endog_prior'
                    'number of active inequalities'
                    'log-MDD(Laplace)'
                    'log-MDD(MHM)'
                    'estimation sample'
                    'number of observations '
                    'number of parameters '
                    'estimation algorithm '
                    'solution algorithm '
                    'start time:'
                    'end time :'
                    'total time:'
                    };
                for icu=1:numel(model_objects)
                    this_ic={model_objects(icu).log_post,model_objects(icu).log_lik,model_objects(icu).log_prior,...
                        model_objects(icu).log_endog_prior,model_objects(icu).numberOfActiveInequalities,...
                        model_objects(icu).log_mdd_laplace,model_objects(icu).log_mdd_mhm,...
                        [model_objects(icu).options.estim_start_date,':',model_objects(icu).options.estim_end_date],...
                        numel(model_objects(icu).varobs(1).value),numel(model_objects(icu).estimated_parameters),...
                        model_objects(icu).options.optimizer,...
                        model_objects(icu).options.solver,nan,nan,nan}';
                    if isfield(model_objects(icu).options,'estim_end_time') && ~isempty(model_objects(icu).options.estim_end_time)
                        t2=model_objects(icu).options.estim_end_time;
                        t1=model_objects(icu).options.estim_start_time;
                        estimation_time=etime(t2,t1);
                        hrs=floor(estimation_time/3600);
                        secs=estimation_time-hrs*3600;
                        mins=floor(secs/60);
                        secs=secs-mins*60;
                        this_ic(end-(2:-1:0))={datestr(t1),datestr(t2),[int2str(hrs),':',int2str(mins),':',int2str(secs)]};
                    end
                    thisTable=[thisTable,[model_names(icu);this_ic]];
                end
                
                thisTable=struct('longtable',false,'table',{thisTable},...
                    'caption','Estimation summary statistics',...
                    'title','Extended estimation output');
                table(obj,thisTable)
            end
        end
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
                myscript=[myscript;{['\section{',reprocess(new_table.title),'}']}];
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
%                 tablestyle=['long',tablestyle];
            end
            myscript=[myscript;{
                ['\begin{',tablestyle,'}[h] \centering']
                ['\begin{tabular}{',repmat('r',1,ncols),'}']
                '\hline\hline'
                theHeader
                '\hline'
                }];
            myscript=[myscript;AllBatch];
            myscript=[myscript;{'\hline'
                '\hline'
                '\end{tabular}'}];
            if ~isempty(new_table.caption)
                myscript=[myscript;{['\caption{',new_table.caption,'}']}];
            end
            myscript=[myscript;{['\end{',tablestyle,'}']}];
            record(obj,myscript)
            obj.clearpage();
            thisobj=obj;
        end
        function thisobj=figure(obj,figure_struct)
            default_figure=struct('name','',...
                'title','',...
                'caption','','angle',0);
            if isempty(obj)
                thisobj=default_figure;
                return
            end
            new_figure=mysetfield(default_figure,figure_struct); %
            if ~isempty(new_figure.name)
                tmpname=new_figure.name;
                angle=new_figure.angle;
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
                if any(isspace(myfigname))
                    myfigname=['"',myfigname,'"'];
                end
                newscript={ %'\newpage'
                    '\begin{tabular}[t]{@{\hspace*{-3pt}}c@{}}'
                    ['\multicolumn{1}{c}{\large\bfseries ',myfigname,'}\\']
                    ['\raisebox{10pt}{\includegraphics[scale=0.9,angle=',...
                    num2str(angle),']{',tmpname,'}}']
                    '\end{tabular}'};
                record(obj,newscript)
                obj.clearpage();
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
            
            retcode=system(['pdflatex ',obj.report_name]);
            if ~retcode
                % run again in order to make sure the references are shown
                system(['pdflatex ',obj.report_name]);
                system(['pdflatex ',obj.report_name]);
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
        function include(obj,filename)
            record(obj,['\input{',filename,'}']);
        end
        function verbatim(obj,varargin)
            record(obj,['\begin{verbatim}';
                varargin(:);
                '\end{verbatim}'])
        end
        function text(obj,CellItems)
            record(obj,CellItems(:))
        end
        function paragraph(obj,varargin)
            text(obj,[' ';varargin(:)])
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
        function enumerate(obj,CellItems)
            for ii=1:numel(CellItems)
                CellItems{ii}=['\item ',CellItems{ii}];
            end
            CellItems=CellItems(:);
            record(obj,['\begin{enumerate}'
                CellItems
                '\end{enumerate}'])
        end
        function itemize(obj,CellItems)
            for ii=1:numel(CellItems)
                CellItems{ii}=['\item ',CellItems{ii}];
            end
            CellItems=CellItems(:);
            record(obj,['\begin{itemize}'
                CellItems
                '\end{itemize}'])
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

function xin=reprocess(xin,precision)
if nargin<2
    precision=4;
end
if isnumeric(xin)
    xin=num2str(xin,precision);
elseif ischar(xin)
    if any(xin=='$')
        return
    end
    xin=strrep(xin,'\_','LouisPergaud');
    xin=strrep(xin,'_','\_');
    xin=strrep(xin,'LouisPergaud','\_');
end
end