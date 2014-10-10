classdef titlepage < rise_report.generic_report
    % titlepage report title page object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.titlepage/addlistener)
    % - [best_title](rise_report.titlepage/best_title)
    % - [delete](rise_report.titlepage/delete)
    % - [eq](rise_report.titlepage/eq)
    % - [findobj](rise_report.titlepage/findobj)
    % - [findprop](rise_report.titlepage/findprop)
    % - [ge](rise_report.titlepage/ge)
    % - [gt](rise_report.titlepage/gt)
    % - [isvalid](rise_report.titlepage/isvalid)
    % - [le](rise_report.titlepage/le)
    % - [lt](rise_report.titlepage/lt)
    % - [ne](rise_report.titlepage/ne)
    % - [notify](rise_report.titlepage/notify)
    % - [reprocess](rise_report.titlepage/reprocess)
    % - [titlepage](rise_report.titlepage/titlepage)
    % - [write](rise_report.titlepage/write)
    %
    % properties
    % -----------
    %
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
        title=''
        date=''
        author=''
        address=''
        email=''
        abstract={}
        latex_date_format=false;
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=titlepage(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            obj.author=cellform(obj.author);
            obj.address=cellform(obj.address);
            obj.email=cellform(obj.email);
            obj.abstract=cellform(obj.abstract);
            b={''};
            if ~isempty(obj.title)
                titel=obj.title; % titel=rise_report.generic_report.best_title(obj.title);
                b={['\title{',titel,'}']};
                the_date=obj.date;
                if isempty(the_date)
                    if obj.latex_date_format
                        the_date='\today';
                        b=[b;the_date];
                    else
                        the_date=datestr(now);
                        b=[b;['\date{',the_date,'}']];
                    end
                else
                    b=[b;['\date{',the_date,'}']];
                end
                n_authors=numel(obj.author);
                if n_authors
                    AfterAuthor=' \\';
                    flag_add=numel(obj.address)==n_authors;
                    flag_email=numel(obj.email)==n_authors;
                    b=[b;'\author{'];
                    for iaut=1:n_authors
                        thisAuthor=obj.author{iaut};
                        if flag_add
                            thisAuthor=[thisAuthor,AfterAuthor];
                            b=[b;thisAuthor];
                            thisAddress=obj.address{iaut};
                            if flag_email
                                thisAddress=[thisAddress,AfterAuthor];
                                b=[b;thisAddress];
                                thisEmail=['\texttt{',obj.email{iaut},'}'];
                                b=[b;thisEmail];
                            else
                                b=[b;thisAddress];
                            end
                        else
                            if flag_email
                                thisAuthor=[thisAuthor,AfterAuthor]; %#ok<*AGROW>
                                b=[b;thisAuthor];
                                thisEmail=['\texttt{',obj.email{iaut},'}'];
                                b=[b;thisEmail];
                            else
                                b=[b;thisAuthor];
                            end
                        end
                        if iaut<n_authors
                            b=[b;'\and'];
                        end
                    end
                    b=[b;'}'];
                end
                % will the following work for the headers?
                b=[b;'\pagestyle{myheadings}'];
                b=[b;['\markright{',titel,'\hfill ',the_date,'\hfill}']];
                b=[b;'\thispagestyle{empty}'];
                b=[b;'\maketitle'];
                if ~isempty(obj.abstract)
                    b=[b;'\begin{abstract}'];
                    for irow=1:numel(obj.abstract)
                        b=[b;obj.abstract{irow}];
                    end
                    b=[b;'\end{abstract}'];
                end
                b=[b;'\pagebreak'];
            end
        end
    end
end

function c=cellform(c)
if ~isempty(c) && ischar(c)
    c=cellstr(c);
end
end