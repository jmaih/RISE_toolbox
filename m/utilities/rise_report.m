function retcode=rise_report(this,report_name,instructions,graphicsPath)

if nargin<4
    graphicsPath='';
end
if ~isempty(graphicsPath)
    graphicsPath=strrep(graphicsPath,'\','/');
    if ~strcmp(graphicsPath(end),'/')
        graphicsPath=[graphicsPath,'/'];
    end
end
default_name='report';
if isempty(report_name)
    report_name=default_name;
end
rise_pdflatex=getappdata(0,'rise_pdflatex');

if ~rise_pdflatex
    error([mfilename,':: cannot generate a report, MIKTEX was not found earlier'])
end
%% create a temporary directory to temporarily save handles
tmpdir='';
tmpdir_flag=false;

thedot=strfind(report_name,'.');
if ~isempty(thedot)
    report_name=report_name(1:thedot-1);
end
if exist([report_name,'.tex'],'file')
    delete([report_name,'.tex'])
end

fid=fopen([report_name,'.tex'],'w');

add_preamble();

% locate the title and add it at the very beginning
default_title={'title';add_title()};
title_loc=find(strcmp(instructions(1,:),'title'), 1);
if ~isempty(title_loc)
    default_title=instructions(:,title_loc);
    instructions=[instructions(:,1:title_loc-1),instructions(:,title_loc+1:end)];
end
instructions=[default_title,instructions];

for iinstr=1:size(instructions,2)
    type=instructions{1,iinstr};
    switch type
        case 'title'
            add_title(instructions{2,iinstr});
        case 'paragraph'
            add_paragraph(instructions{2,iinstr});
        case 'table'
            add_table(instructions{2,iinstr});
        case 'figure'
            add_figure(instructions{2,iinstr});
        case 'equations'
            add_model_equations(instructions{2,iinstr})
        case 'estimation_statistics'
            if isempty(instructions{2,iinstr})
                instructions(:,iinstr)=further_estimation2table();
            end
            add_table(instructions{2,iinstr});
        case 'estimation'
            if isempty(instructions{2,iinstr})
                instructions(:,iinstr)=estimation2table();
                add_table(instructions{2,iinstr});
                fprintf(fid,'%s \n','\newpage');
                instructions(:,iinstr)=further_estimation2table();
            end
            add_table(instructions{2,iinstr});
        otherwise
            error([mfilename,':: unknown type ',type])
    end
    fprintf(fid,'%s \n','\newpage');
end

add_finishing();

fclose(fid);

if exist([report_name,'.pdf'],'file')
    delete([report_name,'.pdf'])
end

retcode=system(['pdflatex ',report_name]);
useless_extensions={'.log','.bbl','.blg','.aux','.*.bak'};
for iext=1:numel(useless_extensions)
    file2delete=[report_name,useless_extensions{iext}];
    if exist(file2delete,'file')
        delete(file2delete)
    end
end

if tmpdir_flag
rmdir(tmpdir,'s'); % use s to remove the contents as well otherwise the operation will not be successful
end
    function add_table(options)
        default_table=struct('title','no title',...
            'table',[],'caption','no caption');
        new_table=mysetfield(default_table,options);
        if isempty(new_table.table)
            return
        end
        fprintf(fid,'%s \n',['\section*{',reprocess(new_table.title),'}']);
        [nrows,ncols]=size(new_table.table);
        AllBatch=cell(nrows,1);
        for irow=1:nrows
            str='';
            for icol=1:ncols
                tmp=reprocess(new_table.table{irow,icol});
% % % % %                 if irow==1
% % % % %                     tmp=['\textbf{',tmp,'}'];
% % % % %                 end
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
        new_nrows=nrows-1;
        max_rows=15;
        number_of_tables=ceil(new_nrows/max_rows);
        additional_string='';
        for itab=1:number_of_tables
            if itab>1
                fprintf(fid,'%s ','Table continued ...');
                fprintf(fid,'%s ','\clearpage');
            end
            if number_of_tables>1
                additional_string=['(',int2str(itab),')'];
            end
            % with this I can cut the table at any time
            fprintf(fid,'%s \n','\begin{table}[h] \centering');
            fprintf(fid,'%s \n',['\begin{tabular}{',repmat('r',1,ncols),'}']);
            fprintf(fid,'%s \n','\hline\hline');
            fprintf(fid,'%s \n',theHeader);
            fprintf(fid,'%s \n','\hline');
            for irow=(itab-1)*max_rows+1:min(new_nrows,itab*max_rows)
                fprintf(fid,'%s \n',AllBatch{irow});
            end
            fprintf(fid,'%s \n','\hline');
            if itab==number_of_tables
                fprintf(fid,'%s \n','\hline');
            end
            fprintf(fid,'%s \n','\end{tabular}');
            fprintf(fid,'%s \n',['\caption{',new_table.caption,additional_string,'}']);
            fprintf(fid,'%s \n','\end{table}');
        end
    end

    function add_preamble()
        preamble_instructions={
            '\documentclass[12pt,landscape]{article}'
            '\usepackage{amsfonts,amsmath,color,graphicx}'%,fullpage
            '\usepackage{geometry}'%[centering]
% %             '\usepackage[space]{grffile}'
            '\begin{document}'
            '\pagestyle{myheadings}'
            };
        if ~isempty(graphicsPath)
            preamble_instructions=[preamble_instructions
            ['\graphicspath{{',graphicsPath,'}}']];
        end
        for jj=1:numel(preamble_instructions)
            fprintf(fid,'%s \n',preamble_instructions{jj});
        end
        fprintf(fid,'%s \n\n','% end of preamble');
    end

    function add_finishing()
        fprintf(fid,'\n%s \n','\end{document}');
    end

    function default_title=add_title(options)
        default_title=struct('title','Simple RISE Report',...
            'address','',...
            'date','',...
            'author','RISE Toolbox',...
            'email','junior.maih@gmail.com');
        if nargin==0
            return
        end
        AfterAuthor=' \\';
        new_title=mysetfield(default_title,options);
        fprintf(fid,'%s \n',['\title{',reprocess(new_title.title),'}']);
        authors=regexp(new_title.author,'&','split');
        address=regexp(new_title.address,'&','split');
        email=regexp(new_title.email,'&','split');
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
        fprintf(fid,'%s \n','\maketitle');
    end

    function add_paragraph(paragraph)
        default_paragraph=struct('text','',...
            'title','no title','itemize',false);
        new_paragraph=mysetfield(default_paragraph,paragraph);
        new_paragraph.text=char(new_paragraph.text);
        itemize=new_paragraph.itemize;
        fprintf(fid,'%s \n',['\section*{',reprocess(new_paragraph.title),'}']);
        string='';
        if itemize
            fprintf(fid,'%s \n','\begin{itemize}');
            string='\item ';
        end
        for ii=1:size(new_paragraph.text,1)
            fprintf(fid,'%s \n',[string,'\large{',new_paragraph.text(ii,:),'}']);
        end
        if itemize
            fprintf(fid,'%s \n','\end{itemize}');
        end
    end

    function add_model_equations(equations)
        default_equations=struct('title','Model equations',...
            'equations',{{this(1).equations.dynamic}});
        new_equations=mysetfield(default_equations,equations);
        if ischar(new_equations.equations)
            new_equations.equations=cellstr(new_equations.equations);
        end
        fprintf(fid,'%s \n',['\section*{',reprocess(new_equations.title),'}']);
        fprintf(fid,'%s \n',' \begin{verbatim}');	%\color{lightgray}
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
                fprintf(fid,'%s \n',myeq(isubeq,:));
            end
            fprintf(fid,'%s \n','');
        end
        fprintf(fid,'%s \n','\end{verbatim}\color{black} ');
    end

    function add_figure(graph)
        default_figure=struct('name','',...
            'title','no title',...
            'caption','no caption','angle',0);
        new_figure=mysetfield(default_figure,graph); %
        if ~isempty(new_figure.name)
            tmpname=new_figure.name;
            angle=new_figure.angle;
            if ishandle(tmpname)
                if ~tmpdir_flag
                    tmpdir=tempname(pwd);
                    mkdir(tmpdir);
                    tmpdir_flag=true;
                end
                new_handle=tmpname;
                tmpname=tempname(tmpdir); 
               angle=rise_saveaspdf(new_handle,tmpname);
            end
            tmpname=strrep(tmpname,'\','/');
            tmpname=strrep(tmpname,'.pdf','');%,'.pdf'
             % quotes and spaces may still not work if the extension .pdf is
            % added
            fprintf(fid,'%s \n','\newpage');
            fprintf(fid,'%s \n','\begin{tabular}[t]{@{\hspace*{-3pt}}c@{}}');
            myfigname=reprocess(new_figure.title);
            if any(isspace(myfigname))
                myfigname=['"',myfigname,'"'];
            end
            fprintf(fid,'%s \n',['\multicolumn{1}{c}{\large\bfseries ',myfigname,'}\\']);
            fprintf(fid,'%s \n',['\raisebox{10pt}{\includegraphics[scale=0.9,angle=',...
                num2str(angle),']{',tmpname,'}}']);
            fprintf(fid,'%s \n','\end{tabular}');%
        end
    end

    function mytable=estimation2table()
        ncases=numel(this);
        type_name='tex_name';
        parnames= {this(1).estimated_parameters.(type_name)}';
        for ic=2:ncases
            if ~isequal(parnames,{this(ic).estimated_parameters.(type_name)}')
                error([mfilename,':: all models should have the same (estimated) parameters'])
            end
        end
        estimated_parameters=[this.estimated_parameters];
        npar=numel(parnames);
        mytable=cell(npar+1,5);
        for ipar=1:npar
            if any(parnames{ipar}=='\')
                if ~strcmp(parnames{ipar}(1),'$')
                    parnames{ipar}=['$',parnames{ipar}];
                end
                if ~strcmp(parnames{ipar}(end),'$')
                    parnames{ipar}=[parnames{ipar},'$'];
                end
            end
        end
        mytable(2:end,1)=parnames;mytable{1,1}='Param Names';
        mytable(2:end,2)={estimated_parameters(:,1).distribution}';mytable{1,2}='Prior distr';
        mytable(2:end,3)=num2cell(vertcat(estimated_parameters(:,1).interval_probability));mytable{1,3}='Prior prob';
        mytable(2:end,4)=num2cell(vertcat(estimated_parameters(:,1).plb));mytable{1,4}='low';
        mytable(2:end,5)=num2cell(vertcat(estimated_parameters(:,1).pub));mytable{1,5}='high';
        for ic=1:ncases
            mytable=[mytable,[{['mode(',int2str(ic),')']};num2cell(vertcat(estimated_parameters(:,ic).mode))]];
        end
        if ncases==1
            mytable=[mytable,[{'mode\_std'};num2cell(vertcat(estimated_parameters(:,1).mode_std))]];
        end
        mytable={'table'
            struct('title','Estimation Results',...
            'table',{mytable},'caption','Prior and Posterior mode')};
    end
    function mytable=further_estimation2table()
        mytable={
            'log-post:'
            'log-lik:'
            'log-prior:'
            'log-endog_prior'
            'numberOfActiveInequalities'
            'log-MDD(Laplace)'
            'log-MDD(MHM)'
            'estimation sample'
            'number of observations '
            'estimation algorithm '
            'solution algorithm '
            'start time:'
            'end time :'
            'total time:'
            };
        for ic=1:numel(this)
            this_ic={this(ic).log_post,this(ic).log_lik,this(ic).log_prior,...
                this(ic).log_endog_prior,this(ic).numberOfActiveInequalities,...
                this(ic).log_mdd_laplace,this(ic).log_mdd_mhm,...
                [this(ic).options.estim_start_date,':',this(ic).options.estim_end_date],...
                numel(this(ic).varobs(1).value),this(ic).options.optimizer,...
                this(ic).options.solver,nan,nan,nan}';
            if isfield(this(ic).options,'estim_end_time') && ~isempty(this(ic).options.estim_end_time)
                t2=this(ic).options.estim_end_time;
                t1=this(ic).options.estim_start_time;
                estimation_time=etime(t2,t1);
                hrs=floor(estimation_time/3600);
                secs=estimation_time-hrs*3600;
                mins=floor(secs/60);
                secs=secs-mins*60;
                this_ic(end-(2:-1:0))={datestr(t1),datestr(t2),[int2str(hrs),':',int2str(mins),':',int2str(secs)]};
            end
            mytable=[mytable,this_ic];
        end
        
        mytable={'table'
            struct('title','Extended estimation output',...
            'table',{mytable},'caption','Estimation summary statistics')};
    end
end

function xin=reprocess(xin,precision)
if nargin<2
    precision=4;
end
if isnumeric(xin)
    xin=num2str(xin,precision);
elseif ischar(xin)
    xin=strrep(xin,'\_','LouisPergaud');
    xin=strrep(xin,'_','\_');
    xin=strrep(xin,'LouisPergaud','\_');
end
end
