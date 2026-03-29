function publish(obj,savefile,varargin)
% PUBLISH - Save the report to a file.
%
% Syntax:
%
%   obj.publish()
%
%   obj.publish(savefile)
%
%   obj.publish(savefile, optName1, optValue1, ...)
%
% Description:
%
%   This function saves the report to a file in either HTML or PDF format.
%   By default, the report is saved as 'report.html'. The user can specify
%   the file name and additional options using the input arguments.
%
% Inputs:
%
%   - obj : The report object.
%
%   - savefile : Optional. Name of the file to save the report (default:
%     'report.html'). The name can be given without extension, in which
%     case the default is html. Only two extensions are allowed : html and
%     pdf.
%
%   - varargin : Optional. Additional options to customize the report
%     generation. 
%
% Options:
%
%   - 'clean_up' : Optional. Option to delete useless files (default: 'all').
%                  Valid values are 'all', 'none', or 'all_but_tex'. This
%                  option works only for pdf files. The 'all_but_tex'
%                  variant avoids the deletion of the tex file through
%                  which the pdf file is generated
%
%   - 'graphics_path' : Optional. Path to the graphics folder (default: empty).
%
% Example:
%
%   report = rprt('Sample Report', 'John Doe', 'portrait');
%   % Add elements to the report...
%   report.publish();                 % Saves the report to 'report.html'
%   report.publish('myreport.pdf');   % Saves the report to 'myreport.pdf'
%
% Notes:
%
%   - If a 'graphics_path' is not provided, the object will create a folder
%     based on the name of the 'savefile'.
%
%   - The saving format is determined by the extension of the 'savefile' name.
%     If no extension is provided, the default format is HTML.

options = update_options();

% Check the number of input arguments

if nargin < 2||isempty(savefile)
    
    savefile='report.html';
    
end

% Extract the rootname and extension from the savefile
[~, rootname, ext] = fileparts(savefile);

% Set default extension to 'html' if not provided
if isempty(ext)
    
    ext = '.html'; % Default extension is 'html' if not provided
    
    savefile = [rootname, ext];
    
end

% provision for pdflatex
swapext=@(x)strrep(x,'.pdf','.tex');

% Open the file for writing
fid = fopen(swapext(savefile), 'w');

if fid == -1
    
    error('Unable to open file for writing: %s', savefile);
    
end

% Determine the format based on the file extension
format = ext;

folders_created={};

if isempty(options.graphics_path)
    
    options.graphics_path=[rootname,'_figures'];%tempname(pwd);
    
    folders_created=[folders_created,options.graphics_path];
    
end

if ~isdir(options.graphics_path) %#ok<ISDIR>
    
    mkdir(options.graphics_path);
    
end

mv=ver; mv=mv(strcmpi({mv.Name},'matlab'));

mv=[mv.Version,mv.Release];

rv=rise.version();

close_file=true;

% Generate the report based on the format
if strcmpi(format, '.html')
    
    frmt = 'html';
    
    toHtml()
    
elseif strcmpi(format, '.pdf')
    
    frmt = 'pdf';
    
    tex_file=toPdf();
    
    pdflatex(tex_file,options,folders_created);
    
    % Close the file
    close_file= ~options.clean_up;
    
else
    
    error('Invalid file format. Only HTML and PDF formats are supported.');
    
end

% Close the file
if close_file
    
    fclose(fid);
    
end

disp(['Report saved as: ', savefile]);


    function tex_file=toPdf()
                
        fprintf(fid,'%% Report Object\n');
        
        fprintf(fid, '\\documentclass[%s,%s,%s]{%s}\n',...
            obj.papersize,obj.points,obj.orientation,obj.documentclass);
        
        % blocked packages
        %-------------------
        blocked_list={'inputenc','fontenc','geometry'};
        fprintf(fid, '\\usepackage[utf8x]{inputenc}\n');% fprintf(fid, '\\usepackage[utf8]{inputenc}\n');
        fprintf(fid, '\\usepackage[T1]{fontenc}\n');
        fprintf(fid, '\\usepackage[nice]{nicefrac}\n');
        fprintf(fid, '\\usepackage[pscoord]{eso-pic}\n');
        fprintf(fid, '\\usepackage[margin=1in]{geometry}\n');
        other_packages={'multirow','adjustbox','natbib','etex','xcolor',...
            'fourier-orns','ccicons','amssymb,amstext,amsbsy,amsopn,amscd,amsxtra,amsthm',...
            'float','color,soul','mathrsfs','bm','lastpage','xstring',...
            'setspace','hyperref','ragged2e','listings',...
            'algorithms/algorithm','algorithms/algorithmic','tikz,pgfplots',...
            'amsmath,verbatim','amsfonts','graphicx','lipsum'};

        packages=union(obj.packages,other_packages);
        % optional packages
        %-------------------
        for ii=1:numel(packages)
            
            packi=packages{ii};
            
            if ~ismember(packi,blocked_list)
                
                fprintf(fid, '\\usepackage{%s}\n',packi);
                
            end
            
        end
        fprintf(fid,'\\graphicspath{{%s}}\n',options.graphics_path);
        fprintf(fid, '\\numberwithin{equation}{section}\n');
        
        fprintf(fid, '\\lstset{\n');
        fprintf(fid, 'basicstyle=\\ttfamily,\n');
        fprintf(fid, 'keywordstyle=\\color{blue},\n');
        fprintf(fid, 'commentstyle=\\color{green},\n');
        fprintf(fid, 'stringstyle=\\color{red},\n');
        fprintf(fid, 'escapechar=|,\n');
        fprintf(fid, '}\n');
        
        fprintf(fid, '\\begin{document}\n');
        
        if isempty(obj.Title)
            
            fprintf(fid,'\\title{}\n');
            
        else
            
            fprintf(fid, '\\title{%s}\n', obj.Title);
            
        end
        
        write_pdf_authors(fid,obj.Author)
        
        fprintf(fid, '\\date{%s\\\\\\bigskip\n', obj.Date);
        fprintf(fid, 'RISE version \\# %s \\& Matlab version %s}\n', rv,mv);
        fprintf(fid, '\\maketitle\n');
        % The possibly some abstract right after maketitle
        
        fprintf(fid,'\\pagebreak\n');
        fprintf(fid,'\\pagestyle{fancy}\n');
        fprintf(fid,'\\fancyhf{}\n');
        
        fprintf(fid, '\\lhead{Generated on %s}\n',obj.Date); % Left header content (leave empty for no content)
        fprintf(fid, '\\rhead{using RISE version %s \\& Matlab version %s}\n',rv,mv); % Header content
        fprintf(fid, '\\lfoot{}\n'); % Left footer content (leave empty for no content)
        fprintf(fid, '\\rfoot{}\n'); % Right footer content (leave empty for no content)
        
        environment='%-------------------------------------------------%';
        
        % Loop over the elements in the Body and write them to the file
        for i = 1:numel(obj.Body)
            
            fprintf(fid,'%s\n',environment);
            
            fprintf(fid,'%% item # %0.0f\n',i);
            
            obj.Body{i}(fid, frmt,options.graphics_path);
            
            fprintf(fid,'%s\n',environment);
            
        end
        
        fprintf(fid, '\\end{document}\n');
        
        tex_file=swapext(savefile);
        
    end


    function toHtml()
        
        fprintf(fid,'<!DOCTYPE html>\n');
        fprintf(fid,'<html>\n');
        fprintf(fid,'<head>\n');
        fprintf(fid,'\t<title>%s</title>\n',obj.Title);
        fprintf(fid,'\t<meta charset="UTF-8">\n');
        fprintf(fid,'\t<script type="text/javascript"\n');
        fprintf(fid,'\t\tsrc="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">\n');
        fprintf(fid,'\t</script>\n');
        if obj.EquationNumber
            fprintf(fid,'\t<style>\n');
            fprintf(fid,'\t.equation {\n');
            fprintf(fid,'\ttext-align: center;\n');
            fprintf(fid,'\tmargin-bottom: 10px;\n');
            fprintf(fid,'\t}\n');
            
            fprintf(fid,'\t.equation-number {\n');
            fprintf(fid,'\tfloat: right;\n');
            fprintf(fid,'\tfont-weight: bold;\n');
            fprintf(fid,'\ttext-decoration: underline;\n');
            fprintf(fid,'\tmargin-left: 10px;\n');
            fprintf(fid,'\t}\n');
            
            fprintf(fid,'\t.equation-expression {\n');
            fprintf(fid,'\ttext-align: center;\n');
            fprintf(fid,'\t}\n');
            fprintf(fid,'\t</style>\n');
            
        end
        fprintf(fid,'</head>\n');
        fprintf(fid,'<body>\n');
        
        % Write the title, author, and date information to the HTML file
        if ~isempty(obj.Title)
            
            fprintf(fid,'\t<h1>%s</h1>\n',obj.Title);
            
        end
        
        write_html_authors(fid, obj.Author)
        
        fprintf(fid,'\t<p>Date: %s</p>\n',obj.Date);
        fprintf(fid,'\t<p>RISE version %s &amp; Matlab version %s</p>\n', rv,mv);
        
        % Loop over the elements in the Body and write them to the file
        for i = 1:numel(obj.Body)
            
            obj.Body{i}(fid, frmt,options.graphics_path);
            
        end
        
        write_footnotes(obj, fid)
        
        % Close the HTML file
        fprintf(fid,'</body>\n');
        fprintf(fid,'</html>\n');
        
        % html needs the generated files to operate and so we do not delete
        % them
        
    end


    function options = update_options()
        % UPDATE_OPTIONS Updates the options structure for the 'figure' method.
        d = {
            'graphics_path', '', @(x) ischar(x), 'graphics_path should be a string'
            
            'clean_up','all',@(x)ismember(x,{'all','none','all_but_tex'}),...
            'clean_up should be "all" or "none" or "all_but_tex"'
            };
        
        options = rprt.set_options(d, mfilename, varargin{:});
    end

end


function write_pdf_authors(fid,authors)

if isempty(authors)
    
    fprintf(fid, '\\author{}\n');
    
    return
    
end

isEmail=all(cellfun(@(x)~isempty(x),authors(:,2)));

na=size(authors,1);

allInstitutes=cell(1,0);

for ia=1:na
    
    instia=authors{ia,3};
    
    if isempty(instia),continue,end
    
    instia=cellstr(instia);
    
    allInstitutes=[allInstitutes,instia]; %#ok<AGROW>
    
    if ia==na
        
        allInstitutes=unique(allInstitutes);
        
    end
    
end

numberOfInstitutions=numel(allInstitutes);

for ia=1:na
    
    fprintf(fid, '\\author%s{%s%s}\n',affil(ia),authors{ia,1},thanks(ia));
    
end

for jj=1:numberOfInstitutions
    
    fprintf(fid, '\\affil[%d]{%s}\n',jj,allInstitutes{jj});
    
end


    function t=thanks(indx)
        
        t='';
        
        if ~isEmail
            
            return
            
        end
        
        t=sprintf('\\thanks{%s}',authors{indx,2});
        
    end


    function af=affil(indx)
        
        af='';
        
        if numberOfInstitutions==0
            
            return
            
        end
        
        insti=cellstr(authors{indx,3});
        
        ni=numel(insti);
        
        tmp=zeros(1,ni);
        
        for ii=1:ni
            
            tmp(ii)=find(strcmp(insti{ii},allInstitutes));
            
        end
        
        tmp=num2cell(sort(tmp));
        
        tmp=cellfun(@int2str,tmp,'uniformOutput',false);
        
        tmp=cell2mat(strcat(tmp,','));
        
        af=['[',tmp(1:end-1),']'];
        
    end

end


function pdflatex(filename,options,folders_created)
% INTERNAL FUNCTION
%

report_name=parser.remove_file_extension(filename);

utils.latex.pdflatex(report_name)

% delete all useless files
useless_extensions={'.log','.bbl','.blg','.aux','.bak','.out'};


if ~strcmp(options.clean_up,'none')
    
    if ~strcmp(options.clean_up,'all_but_tex')
        
        useless_extensions=[useless_extensions,'.tex'];
        
    end
    
    for ifolder=1:numel(folders_created)
        
        rmdir(folders_created{ifolder},'s');
        
    end
    
    for iext=1:numel(useless_extensions)
        
        file2delete=[report_name,useless_extensions{iext}];
        
        if exist(file2delete,'file')
            
            delete(file2delete)
            
        end
        
    end
    
end

end


function write_html_authors(fid, authors)

if isempty(authors)
    
    fprintf(fid, '<author></author>\n');
    
    return
    
end

na = size(authors, 1);

for ia = 1:na
    
    fprintf(fid, '<author>\n');
    
    fprintf(fid, '    <name>%s</name><br>\n', authors{ia, 1});
    
    fprintf(fid, '    <email>%s</email><br>\n', authors{ia, 2});
    
    fprintf(fid, '    <affiliations>\n');
    
    write_affiliations(fid, authors{ia, 3});
    
    fprintf(fid, '    </affiliations><br>\n');
    
    fprintf(fid, '</author><br><br>\n');
    
end

    function write_affiliations(fid, affiliations)
        
        if ischar(affiliations)
            
            affiliations=cellstr(affiliations);
            
        end
        
        for iaff=1:numel(affiliations)
            
            fprintf(fid, '        <affiliation>%s</affiliation><br>\n', affiliations{iaff});
            
        end
        
    end

end