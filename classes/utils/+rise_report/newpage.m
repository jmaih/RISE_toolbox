classdef newpage < rise_report.generic_report
    % newpage report new page object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.newpage/addlistener)
    % - [best_title](rise_report.newpage/best_title)
    % - [delete](rise_report.newpage/delete)
    % - [eq](rise_report.newpage/eq)
    % - [findobj](rise_report.newpage/findobj)
    % - [findprop](rise_report.newpage/findprop)
    % - [ge](rise_report.newpage/ge)
    % - [gt](rise_report.newpage/gt)
    % - [isvalid](rise_report.newpage/isvalid)
    % - [le](rise_report.newpage/le)
    % - [lt](rise_report.newpage/lt)
    % - [ne](rise_report.newpage/ne)
    % - [newpage](rise_report.newpage/newpage)
    % - [notify](rise_report.newpage/notify)
    % - [reprocess](rise_report.newpage/reprocess)
    % - [write](rise_report.newpage/write)
    %
    % properties
    % -----------
    %
    % - [thispagestyle] -
    % - [batch] -
    % - [id] -
    properties(Constant)
        % plain prints the page numbers on the bottom of the page, in the
        % middle of the footer. This is the default page style.
        % headings prints the current chapter heading and the page number
        % in the header on each page, while the footer remains empty. (This
        % is the style used in this document)
        % empty sets both the header and the footer to be empty.
        thispagestyle=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=newpage(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            if nargin
                if ~any(strcmp(varargin{2},rise_report.page_styles()))
                    error([varargin{2},' is not a valid page style'])
                end
            end
        end
        function b = get.batch(obj)
            b={' ';'\newpage'};
            if ~isempty(obj.thispagestyle)
                b=[
                    b
                    ['\thispagestyle{',obj.thispagestyle,'}']
                    ];
            end
        end
    end
end
