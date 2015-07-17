classdef report < rise_report.titlepage
    % report Reporting system
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.report/addlistener)
    % - [best_title](rise_report.report/best_title)
    % - [chapter](rise_report.report/chapter)
    % - [cleardoublepage](rise_report.report/cleardoublepage)
    % - [clearpage](rise_report.report/clearpage)
    % - [delete](rise_report.report/delete)
    % - [enumerate](rise_report.report/enumerate)
    % - [eq](rise_report.report/eq)
    % - [figure](rise_report.report/figure)
    % - [findobj](rise_report.report/findobj)
    % - [findprop](rise_report.report/findprop)
    % - [footnote](rise_report.report/footnote)
    % - [ge](rise_report.report/ge)
    % - [gt](rise_report.report/gt)
    % - [include](rise_report.report/include)
    % - [isvalid](rise_report.report/isvalid)
    % - [itemize](rise_report.report/itemize)
    % - [le](rise_report.report/le)
    % - [lt](rise_report.report/lt)
    % - [ne](rise_report.report/ne)
    % - [newpage](rise_report.report/newpage)
    % - [notify](rise_report.report/notify)
    % - [pagebreak](rise_report.report/pagebreak)
    % - [paragraph](rise_report.report/paragraph)
    % - [publish](rise_report.report/publish)
    % - [quotation](rise_report.report/quotation)
    % - [report](rise_report.report/report)
    % - [reprocess](rise_report.report/reprocess)
    % - [section](rise_report.report/section)
    % - [subparagraph](rise_report.report/subparagraph)
    % - [subsection](rise_report.report/subsection)
    % - [subsubsection](rise_report.report/subsubsection)
    % - [table](rise_report.report/table)
    % - [text](rise_report.report/text)
    % - [verbatim](rise_report.report/verbatim)
    % - [write](rise_report.report/write)
    %
    % properties
    % -----------
    %
    % - [name] -
    % - [graphicspath] -
    % - [documentclass] -
    % - [packages] -
    % - [orientation] -
    % - [points] -
    % - [papersize] -
    % - [fleqn] -
    % - [leqno] -
    % - [titlepage] -
    % - [onecolumn] -
    % - [twoside] -
    % - [openright] -
    % - [title] -
    % - [date] -
    % - [author] -
    % - [address] -
    % - [email] -
    % - [abstract] -
    % - [latex_date_format] -
    % - [batch] -
    % - [id] -
    properties
        name = 'rise_report_default_name'
        graphicspath='';
        documentclass='article'
        packages={'graphicx','amsmath','geometry','amsfonts',...
            'color','hyperref','longtable','float'}
        %----------------------------------------
        orientation='landscape'
        points='12pt'
        papersize='letterpaper'
        % Typesets displayed formulae left-aligned instead of centred.
        fleqn = ''
        % numbering of formulae on the lhs instead of rhs
        leqno = ''
        % titlepage, notitlepage Specifies whether a new page should be
        % started after the document title or not. The article class does
        % not start a new page by default, while report and book do.
        titlepage = ''
        % onecolumn, twocolumn Instructs LATEX to typeset the document in
        % one column or two columns.
        onecolumn = ''
        % twoside, oneside Specifies whether double or single sided output
        % should be generated.
        twoside = ''
        % openright, openany Makes chapters begin either only on right hand
        % pages or on the next page available.
        openright = ''
    end
    properties(Hidden)
        log={}
        nitems=0
        iter=0
        ncalls=0
        folders_created={}
        table_count=0
        figure_count=0
    end
    methods
        function obj=report(varargin)
            n=length(varargin);
            mypackages={'graphicx','amsmath','geometry','amsfonts',...
                'color','hyperref','longtable','float'};
            own_props={
                'name',@(x)ischar(x)
                'documentclass',@(x)ismember(x,{'article','report','book',...
                'proc','minimal','slides'})
                'orientation',@(x)ismember(x,{'portrait','landscape'})
                'points',@(x)ismember(x,{'10pt','11pt','12pt'})
                'papersize',@(x)ismember(x,{'a4paper','letterpaper',...
                'a5paper','b5paper','executivepaper','legalpaper'})
                'packages',@(x)all(cellfun(@(x)ismember(x,mypackages),'uniformOutput',false))
                'fleqn',@(x)strcmp(x,'fleqn')
                'leqno',@(x)strcmp(x,'leqno')
                'titlepage',@(x)any(strcmp(x,{'titlepage','notitlepage'}))
                'onecolumn',@(x)any(strcmp(x,{'onecolumn','twocolumn'}))
                'twoside',@(x)any(strcmp(x,{'twoside','oneside'}))
                'openright',@(x)any(strcmp(x,{'openright','openany'}))
                'graphicspath',@(x)ischar(x) && isdir(x)
                };
            nop=size(own_props,1);
            own_vals=cell(nop,1);
            discard=false(1,n);
            for op=1:nop
                id=find(strcmp(own_props{op,1},varargin));
                if ~isempty(id)
                    discard([id,id+1])=true;
                    own_vals{op}=varargin{id+1};
                end
            end
            varargin=varargin(~discard);
            obj=obj@rise_report.titlepage(varargin{:});
            for op=1:nop
                if ~isempty(own_vals{op})
                    assert(own_props{op,2}(own_vals{op}))
                    obj.(own_props{op,1})=own_vals{op};
                end
            end
            obj.id=0;
            if isempty(obj.graphicspath)
                obj.graphicspath=[obj.name,'_figures'];%tempname(pwd);
                obj.folders_created=[obj.folders_created,obj.graphicspath];
            end
            if ~isdir(obj.graphicspath)
                mkdir(obj.graphicspath);
            end
        end
        function text(obj,varargin)
            record(obj,'text',varargin{:})
        end
        function enumerate(obj,varargin)
            record(obj,'enumerate',varargin{:})
        end
        function itemize(obj,varargin)
            record(obj,'itemize',varargin{:})
        end
        function footnote(obj,varargin)
            record(obj,'footnote',varargin{:})
        end
        function chapter(obj,varargin)
            record(obj,'chapter',varargin{:})
        end
        function subparagraph(obj,varargin)
            record(obj,'subparagraph',varargin{:})
        end
        function paragraph(obj,varargin)
            record(obj,'paragraph',varargin{:})
        end
        function subsubsection(obj,varargin)
            record(obj,'subsubsection',varargin{:})
        end
        function subsection(obj,varargin)
            record(obj,'subsection',varargin{:})
        end
        function section(obj,varargin)
            record(obj,'section',varargin{:})
        end
        function newpage(obj,varargin)
            record(obj,'newpage',varargin{:})
        end
        function include(obj,varargin)
            record(obj,'include',varargin{:})
        end
        function verbatim(obj,varargin)
            record(obj,'verbatim',varargin{:})
        end
        function clearpage(obj)
            record(obj,'clearpage')
        end
        function cleardoublepage(obj)
            record(obj,'cleardoublepage')
        end
        function pagebreak(obj)
            record(obj,'pagebreak')
        end
        function table(obj,varargin)
            record(obj,'table',varargin{:})
        end
        function figure(obj,varargin)
            record(obj,'figure',varargin{:})
        end
        function quotation(obj,varargin)
            record(obj,'quotation',varargin{:})
        end
        %-----------------------------------
        %         function graph(obj,varargin)
        %             record(obj,'graph',varargin{:})
        %         end
        %         function series(obj,varargin)
        %             record(obj,'series',varargin{:})
        %         end
        %         function part(obj,varargin)
        %             record(obj,'part',varargin{:})
        %         end
        %-----------------------------------
        function publish(obj,varargin)
            n=length(varargin);
            if rem(n,2)
                error([mfilename,':: arguments must come in pairs'])
            end
            default_options=struct('write2disk',true,'clean_up',true);
            
            fields=fieldnames(default_options);
            for iarg=1:2:n
                if ~any(strcmp(varargin{iarg},fields))
                    error([varargin{iarg},' is not a valid of option for publish'])
                end
                if ~islogical(varargin{iarg+1})
                    error('Each second argument should be a logical (true or false)')
                end
                default_options.(varargin{iarg})=varargin{iarg+1};
            end
            write2disk=default_options.write2disk;
            clean_up=default_options.clean_up;
            report_name=parser.remove_file_extension(obj.name);
            do_write_up()
            
            if write2disk
                utils.latex.pdflatex(report_name)
            end
            % delete all useless files
            useless_extensions={'.log','.bbl','.blg','.aux','.bak','.out'};
            if clean_up
                useless_extensions=[useless_extensions,'.tex'];
                for ifolder=1:numel(obj.folders_created)
                    rmdir(obj.folders_created{ifolder},'s');
                end
            end
            if write2disk
                for iext=1:numel(useless_extensions)
                    file2delete=[report_name,useless_extensions{iext}];
                    if exist(file2delete,'file')
                        delete(file2delete)
                    end
                end
            end
            
            function do_write_up()
                if write2disk
                    [fid, msg] = fopen([report_name,'.tex'],'w');
                    if fid == -1
                        error([mfilename,':: ',msg]);
                    end
                else
                    fid=1;
                end
                add_preamble();
                add_title_page();
                add_body();
                add_finish();
                %-------------------------------------
                
                function add_finish()
                    fprintf(fid, '\\end{document}\n');
                    fprintf(fid, '%% End Report Object\n');
                    if write2disk
                        status = fclose(fid);
                        if status == -1
                            error([mfilename,':: closing %s\n'], report_name);
                        end
                    end
                    disp('Finished Writing Report!');
                end
                
                function add_body()
                    for iter_=1:obj.iter
                        fprintf(1,'now writing item # %0.0f of %0.0f\n',iter_,obj.iter);
                        write(obj.log{iter_},fid)
                    end
                end
                
                function add_title_page()
                    write(obj,fid)
                end
                
                function add_preamble()
                    fprintf(fid,'%% Report Object\n');
                    doc_attributes=setdiff(properties(rise_report.report),...
                        properties(rise_report.titlepage));
                    doc_attributes=setdiff(doc_attributes,...
                        {'name','documentclass','packages',...
                        'graphicspath','clean_up'});
                    nat=numel(doc_attributes);
                    attrib_vals=cell(1,nat);
                    discard=false(1,nat);
                    for iat=1:nat
                        if isempty(obj.(doc_attributes{iat}))
                            discard(iat)=true;
                        else
                            attrib_vals{iat}=obj.(doc_attributes{iat});
                        end
                    end
                    attrib_vals=attrib_vals(~discard);
                    nat=sum(~discard);
                    tmp=cell2mat(strcat(repmat({'%s'},1,nat),','));
                    fprintf(fid, ['\\documentclass[',tmp(1:end-1),']{%s}\n'],...
                        attrib_vals{:},obj.documentclass);
                    packages_=obj.packages;
                    if ~isempty(packages_)
                        packages_=cell2mat(strcat(packages_(:)',','));
                        fprintf(fid, '\\usepackage{%s}\n',packages_(1:end-1));
                    end
                    fprintf(fid,'\\graphicspath{{%s}}\n',obj.graphicspath);
                    fprintf(fid, '\\numberwithin{equation}{section}\n');
                    %                 fprintf(fid, '\\usepackage[margin=%f%s',...
                    %                     obj.margin, obj.marginUnit);
                    %                 if strcmpi(obj.orientation, 'landscape')
                    %                     fprintf(fid, ',landscape');
                    %                 end
                    %                 fprintf(fid, ']{geometry}\n');
                    %
                    %                     fprintf(fid, '\\usepackage{pdflscape, pgf, booktabs}\n');
                    %                     fprintf(fid, ['\\makeatletter\n' ...
                    %                         '\\def\\blfootnote{\\gdef\\@thefnmark{}\\@footnotetext}\n' ...
                    %                         '\\makeatother\n']);
                    %
                    %                     fprintf(fid, '\\usepackage{pgfplots}\n');
                    %
                    %                     fprintf(fid, '\\usepackage{color, colortbl}\n');
                    %                     fprintf(fid, '\\definecolor{LightCyan}{rgb}{0.88,1,1}\n');
                    %                     fprintf(fid, '\\definecolor{Gray}{gray}{0.9}\n');
                    
                    %                     if obj.showDate
                    %                         fprintf(fid, '\\usepackage{fancyhdr, datetime}\n');
                    %                         fprintf(fid, '\\newdateformat{reportdate}{\\THEDAY\\ \\shortmonthname\\ \\THEYEAR}\n');
                    %                         fprintf(fid, '\\pagestyle{fancy}\n');
                    %                         fprintf(fid, '\\renewcommand{\\headrulewidth}{0pt}\n');
                    %                         fprintf(fid, '\\renewcommand{\\footrulewidth}{0.5pt}\n');
                    %                         fprintf(fid, '\\rfoot{\\scriptsize\\reportdate\\today\\ -- \\currenttime}\n');
                    %                     end
                    
                    % May not need these.....
                    %                     fprintf(fid, '\\renewcommand{\\textfraction}{0.05}\n');
                    %                     fprintf(fid, '\\renewcommand{\\topfraction}{0.8}\n');
                    %                     fprintf(fid, '\\renewcommand{\\bottomfraction}{0.8}\n');
                    %                     fprintf(fid, '\\usepackage[Export,PGF]{adjustbox}\n');
                    %                     fprintf(fid, '\\setlength{\\parindent}{0in}\n');
                    %                     fprintf(fid, '\\newlength\\sectionheight\n');
                    fprintf(fid, '\\begin{document}\n');
                    %                     fprintf(fid, '\\centering\n');
                    %                     fprintf(fid, '\\pagestyle{myheadings}\n');
                end
            end
        end
    end
end

function record(obj,method,varargin)
newitem=rise_report.(method)(varargin{:});
if strcmp(method,'figure')
    obj.figure_count=obj.figure_count+1;
    if isempty(newitem.graphicspath)
        newitem.graphicspath=obj.graphicspath;
    end
    newitem.figure_number=obj.figure_count;
elseif strcmp(method,'table')
    obj.table_count=obj.table_count+1;
    newitem.table_number=obj.table_count;
end
obj.ncalls=obj.ncalls+1;
if ~isempty(newitem)
    obj.iter=obj.iter+1;
    if obj.iter>=obj.nitems
        incmnt=50;
        obj.log=[obj.log;cell(incmnt,1)];
        obj.nitems=obj.nitems+incmnt;
    end
    newitem.id=obj.ncalls;
    obj.log{obj.iter}=newitem;
end
end

