function rise_report(this,report_name,instructions)

default_name='report';
if isempty(report_name)
    report_name=default_name;
end

thedot=strfind(report_name,'.');
if ~isempty(thedot)
    report_name=report_name(1:thedot-1);
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
        case 'estimation'
            if isempty(instructions{2,iinstr})
                instructions(:,iinstr)=estimation2table();
                add_table(instructions{2,iinstr});
                fprintf(fid,'%s \n','\newpage');
                instructions(:,iinstr)=further_estimation2table();
                add_table(instructions{2,iinstr});
            else
                add_table(instructions{2,iinstr});
            end
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
system(['pdflatex ',report_name])
useless_extensions={'.log','.bbl','.blg','.aux','.*.bak'};
for iext=1:numel(useless_extensions)
    file2delete=[report_name,useless_extensions{iext}];
    if exist(file2delete,'file')
        delete(file2delete)
    end
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
                if irow==1
                    str=[str,['\textbf{',tmp,'}']];
                else
                    str=[str,tmp];
                end
                if icol<ncols
                    str=[str,' &'];
                end
            end
            if irow<nrows
                str=[str,' \\ \hline'];
            else
                str=[str,' \\ \hline\hline'];
            end
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
            fprintf(fid,'%s \n',['\begin{tabular}{|',repmat('c|',1,ncols),'}']);
            fprintf(fid,'%s \n','\hline\hline');
            fprintf(fid,'%s \n',theHeader);
            for irow=(itab-1)*max_rows+1:min(new_nrows,itab*max_rows)
                fprintf(fid,'%s \n',AllBatch{irow});
            end
            fprintf(fid,'%s \n','\end{tabular}');
            fprintf(fid,'%s \n',['\caption{',new_table.caption,additional_string,'}']);
            fprintf(fid,'%s \n','\end{table}');
        end
    end

    function add_preamble()
        preamble_instructions={
            '\documentclass[12pt]{article}'
            '\usepackage{amsfonts,amsmath,color,graphicx}'
            '\usepackage[landscape]{geometry}'
            '\begin{document}'
            '\pagestyle{myheadings}'
            };
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
            'author','RISE Toolbox');
        if nargin==0
            return
        end
        AfterAuthor=' \\';
        new_title=mysetfield(default_title,options);
        fprintf(fid,'%s \n',['\title{',reprocess(new_title.title),'}']);
        fprintf(fid,'%s \n',['\author{',new_title.author,AfterAuthor]);
        fprintf(fid,'%s \n',[new_title.address,'}']);
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
            fprintf(fid,'%s \n',[string,new_paragraph.text(ii,:)]);
        end
        if itemize
            fprintf(fid,'%s \n','\end{itemize}');
        end
    end

    function add_model_equations(equations)
        default_equations=struct('title','Model equations',...
            'equations',{{this.equations.dynamic}});
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
            'caption','no caption');
        new_figure=mysetfield(default_figure,graph); %
        if ~isempty(new_figure.name)
            fprintf(fid,'%s \n','\newpage');
%             fprintf(fid,'%s \n','\[');
            fprintf(fid,'%s \n','\begin{tabular}[t]{@{\hspace*{-3pt}}c@{}}');
            fprintf(fid,'%s \n',['\multicolumn{1}{c}{\large\bfseries ',reprocess(new_figure.title),'}\\']);
            % if ishandle create figure and delete after compile. the
            % advantage is that we can portrait or landscape the picture
            % depending on the format of the general report
%             fprintf(fid,'%s \n','\begin{figure}[h]');
%             fprintf(fid,'%s \n','\begin{center}');
            fprintf(fid,'%s \n',['\raisebox{10pt}{\includegraphics[scale=0.85]{',new_figure.name,'}}']);
%             fprintf(fid,'%s \n','\end{center}');
%             fprintf(fid,'%s \n','\end{figure}');
            fprintf(fid,'%s \n','\end{tabular}');%
%             fprintf(fid,'%s \n','\]');
        end
    end

    function mytable=estimation2table()
        estimated_parameters=this.estimated_parameters;
        npar=numel(estimated_parameters);
        mytable=cell(npar+1,7);
        mytable(2:end,1)={estimated_parameters.name}';mytable{1,1}='Param Names';
        mytable(2:end,2)={estimated_parameters.distribution}';mytable{1,2}='Prior distr';
        mytable(2:end,3)=num2cell(vertcat(estimated_parameters.interval_probability));mytable{1,3}='Prior prob';
        mytable(2:end,4)=num2cell(vertcat(estimated_parameters.plb));mytable{1,4}='low';
        mytable(2:end,5)=num2cell(vertcat(estimated_parameters.pub));mytable{1,5}='high';
        mytable(2:end,6)=num2cell(vertcat(estimated_parameters.mode));mytable{1,6}='mode';
        mytable(2:end,7)=num2cell(vertcat(estimated_parameters.mode_std));mytable{1,7}='mode\_std';
        mytable={'table'
            struct('title','Estimation Results',...
            'table',{mytable},'caption','Prior and Posterior mode')};
    end
    function mytable=further_estimation2table()
        mytable={
            'log-post:',this.log_post
            'log-lik:',this.log_lik
            'log-prior:',this.log_prior
            'log-endog_prior',this.log_endog_prior
            'nonlcon_viol_penalty',this.nonlcon_viol_penalty
            'log-MDD(Laplace)', this.log_mdd_laplace
            'log-MDD(MHM)',this.log_mdd_mhm
            'estimation sample', [this.options.estim_start_date,':',this.options.estim_end_date]
            'number of observations ',numel(this.varobs(1).value)
            'estimation algorithm ', this.options.optimizer
            'solution algorithm ', this.options.solver
            };
        
        if isfield(this.options,'estim_end_time') && ~isempty(this.options.estim_end_time)
            t2=this.options.estim_end_time;
            t1=this.options.estim_start_time;
            estimation_time=etime(t2,t1);
            hrs=floor(estimation_time/3600);
            secs=estimation_time-hrs*3600;
            mins=floor(secs/60);
            secs=secs-mins*60;
            mytable=[mytable
                {
                'start time:',datestr(t1)
                'end time :',datestr(t2)
                'total time:',[int2str(hrs),':',int2str(mins),':',int2str(secs)]
                }];
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
